/**
    SPORTA: Spot Detection & Screening of X-Ray Diffraction Images​
    Copyright (c) 2023, Denis SPasyuk, Aaryan Patel

    MIT LICENSE

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    Acknowledgments:
    Canadian Light Source 
    Industrial Science
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <unistd.h>
#include <bshuf_h5filter.h>
#include <math.h>
#include <dirent.h>
#include <regex.h>
#include <ctype.h>
#include <sys/wait.h>

#define UNCLASSIFIED -1
#define UNDEFINED -1
#define NOISE -2
#define CORE_POINT 0

typedef struct {
    double *data;
    int width;
    int height;
    double** pixels;
    double preintensity;
    int* mask; 
} Image;

typedef struct {
    double distance;
    double wavelength;
    double beam_center[2];
    double size[2];
    double psize;
    double threshold;
    int radius_90_percent;
    double resolution;
    int num_spots;
    double intensity;
    int ice;
    double max_distance;
    double snr;
    double postintensity;
    double resolution_95; // New metric
} Hdf5Header;

typedef struct {
    double thresfactor; //intensity scaling factor
    double proximity; //pixel exclusion radius
    int index; // frame number
    bool output; // frame number
    char *filename; //data file name
    char *master; //master file name
    bool traversing; //traversing directories 
    char *data;
    float percent;
    char *searchdir;
    bool gnu;
    double gaussian;
    char path;
    int min_pixel_size;
} config;

typedef struct {
    int label;
    int num_pixels;
    int* x_coords;
    int* y_coords;
    double intensity;
    int cluster_id;
    double reachability_distance;
    double core_distance;
} Spot;


typedef struct PositionSet {
    double startPosition;
    double endPosition;
    struct PositionSet* next;
} PositionSet;

double calculateresolution(double radius);

PositionSet* positionSets = NULL;
Hdf5Header header;
Image img;
config cfg;

// directions for 4-connectivity (up, down, left, right)
int dx[4] = {-1, 1, 0, 0};
int dy[4] = {0, 0, -1, 1};

void printconfig(config cfg) {
    printf("\033[1;33mFilename: %s\033[0m\n", cfg.filename);
    printf("Master: %s\n", cfg.master);
    printf("Proximity: %f\n", cfg.proximity);
    printf("Index: %d\n", cfg.index);
    printf("Output: %d\n", cfg.output);
    printf("Traversing: %s\n", cfg.traversing ? "true" : "false");
    printf("Data: %s\n", cfg.data);
    printf("Search Directory: %s\n", cfg.searchdir);
}

void printheader(const Hdf5Header* header) {
    printf("=================================HEADER============================== \n");
    printf("Detector Distance: %.4f m\n", header->distance);
    printf("Wavelength:        %.4f Å\n", header->wavelength);
    printf("Beam Center:       (%.2f, %.2f) pixels\n", header->beam_center[0], header->beam_center[1]);
    printf("Image Size:        (%.0f, %.0f) pixels\n", header->size[0], header->size[1]);
    printf("Pixel Size:        %.6f m\n", header->psize);
    printf("==================================INFO=============================== \n");
    printf("Avg Raw Intensity:       %.2f\n", img.preintensity);
    printf("Avg Spot Intensity:      %.2f\n", header->postintensity);
    printf("Threshold Value:         %.2f\n", header->threshold);
    printf("Signal Quality (LogSNR): %.2f\n", header->snr);
    
    // Calculate values in mm
    double max_radius_mm = header->max_distance * header->psize * 1000.0;
    printf("Farthest Spot Radius:    %.2f mm (%.2f pixels)\n", max_radius_mm, header->max_distance);
    printf("\033[1;31mBest Spot Resolution:    %.2f Å\033[0m\n", header->resolution);
    if(header->resolution_95 > 0) {
        printf("Resolution (95%%):       %.2f Å\n", header->resolution_95);
    }
    
    // Calculate Edge Resolution
    // Explore all 4 corners to find max radius
    double h = header->size[0];
    double w = header->size[1];
    double bc0 = header->beam_center[0];
    double bc1 = header->beam_center[1];
    
    double d1 = sqrt(pow(0 - bc0, 2) + pow(0 - bc1, 2));
    double d2 = sqrt(pow(h - bc0, 2) + pow(0 - bc1, 2));
    double d3 = sqrt(pow(0 - bc0, 2) + pow(w - bc1, 2));
    double d4 = sqrt(pow(h - bc0, 2) + pow(w - bc1, 2));
    
    double edge_radius = d1;
    if (d2 > edge_radius) edge_radius = d2;
    if (d3 > edge_radius) edge_radius = d3;
    if (d4 > edge_radius) edge_radius = d4;
    
    double edge_res = calculateresolution(edge_radius);
    
    printf("Detector Edge Resolution: %.2f Å\n", edge_res);
    printf("\033[1;31mSpot Count:              %d\033[0m\n", header->num_spots);
    
    PositionSet* currentSet = positionSets;
    int groupNumber = 1;

    if (currentSet == NULL){
        // Check if ice was flagged by hard filter even if no rings were added to list
        if (header->ice) {
             printf("Ice Ring: DETECTED (Hard Filter)\n");
        } else {
             printf("Ice Ring: NONE\n");
        }
    }
    else{
        while (currentSet != NULL) {
            printf("Ice Ring %d: %.2f-%.2f Å\n", groupNumber, currentSet->startPosition, currentSet->endPosition);
            currentSet = currentSet->next;
            groupNumber++;
        }

        while (positionSets != NULL) {
            PositionSet* temp = positionSets;
            positionSets = positionSets->next;
            free(temp);
        }
    }
}

void mainprint(){
    printf("=================================DATASET============================= \n");
    printconfig(cfg);
    printheader(&header);
    printf("===================================================================== \n");
    
}

void noteprint(){
}

void savepgm(const char* filename, double* data, int width, int height) {
    char output_filename[64];
    snprintf(output_filename, sizeof(output_filename), "%s.pgm", filename);
    FILE* file = fopen(output_filename, "w");
    if (file == NULL) {
        printf("Error opening file for writing: %s\n", output_filename);
        return;
    }

    // Write the PGM header
    fprintf(file, "P2\n");
    fprintf(file, "%d %d\n", width, height);
    fprintf(file, "255\n");
    if (strcmp(filename, "pre-output") == 0){
        // Write the inverted image data
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int value = (int)(data[i * width + j]);
                // Clamp to 0-255
                if (value < 0) value = 0;
                if (value > 255) value = 255;
                // Invert
                fprintf(file, "%d ", 255 - value);
            }
            fprintf(file, "\n");
        }
    }
    else {
        // Write the inverted image data for post-output
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int value = (int)(data[i * width + j]);
                // Clamp to 0-255
                if (value < 0) value = 0;
                if (value > 255) value = 255;
                // Invert
                fprintf(file, "%d ", 255 - value);
            }
            fprintf(file, "\n");
        }
    }
    fclose(file);
}


double* decompresshdf5(const char* filename, const char* dataset_path, int index) {
    H5Eset_auto(H5E_DEFAULT, (H5E_auto_t)H5Eprint, stderr);
    herr_t status = bshuf_register_h5filter();
    if (status < 0) {
        fprintf(stderr, "Error registering Bitshuffle filter\n");
        return NULL;
    }

    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return NULL;
    }

    hid_t dataset_id = H5Dopen(file_id, dataset_path, H5P_DEFAULT);
    if (dataset_id < 0) {
        fprintf(stderr, "Error opening dataset: %s\n", dataset_path);
        H5Fclose(file_id);
        return NULL;
    }

    // Check if the dataset is compressed with Bitshuffle
    // hid_t filter_id;
    // unsigned int flags;
    hid_t dataspace_id = H5Dget_space(dataset_id);
    hsize_t dims[3];
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    int num_images = (int)dims[0];

    if (index >= num_images) {
        fprintf(stderr, "Index out of bounds\n");
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
        return NULL;
    }

    hsize_t start[3] = {index, 0, 0};
    hsize_t count[3] = {1, dims[1], dims[2]};
    hid_t mem_dataspace = H5Screate_simple(3, count, NULL);
    double* data_out = (double*)malloc(dims[1] * dims[2] * sizeof(double));

    if (data_out == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        H5Sclose(mem_dataspace);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
        return NULL;
    }

    if (H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
        fprintf(stderr, "Error selecting hyperslab\n");
        H5Sclose(mem_dataspace);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
        free(data_out);
        return NULL;
    }

    if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, mem_dataspace, dataspace_id, H5P_DEFAULT, data_out) < 0) {
        fprintf(stderr, "Error reading data\n");
        H5Sclose(mem_dataspace);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
        free(data_out);
        return NULL;
    }

    H5Sclose(mem_dataspace);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    return data_out;
}

void readhdf5(const char* filename){
    hid_t file_id, detector_id, beam_id, distance_id, detectorX_size,
        detectorY_size, wavelength_id, beamX_id, beamY_id, detectorS_id, psize_id;
    
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    // data_id = H5Gopen(file_id, "/entry/data", H5P_DEFAULT);
    detector_id = H5Gopen(file_id, "/entry/instrument/detector", H5P_DEFAULT);
    detectorS_id = H5Gopen(file_id, "/entry/instrument/detector/detectorSpecific", H5P_DEFAULT);
    beam_id = H5Gopen(file_id, "/entry/instrument/beam", H5P_DEFAULT);
    
    detectorX_size = H5Dopen(detectorS_id, "x_pixels_in_detector", H5P_DEFAULT);
    H5Dread(detectorX_size, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &header.size[0]);  
    H5Dclose(detectorX_size);
    detectorY_size = H5Dopen(detectorS_id, "y_pixels_in_detector", H5P_DEFAULT);
    H5Dread(detectorY_size, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &header.size[1]);  
    H5Dclose(detectorY_size);
    distance_id = H5Dopen2(detector_id, "detector_distance", H5P_DEFAULT);
    H5Dread(distance_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &header.distance);  
    H5Dclose(distance_id);
    wavelength_id = H5Dopen2(beam_id, "incident_wavelength", H5P_DEFAULT);
    H5Dread(wavelength_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &header.wavelength);
    H5Dclose(wavelength_id);
    beamX_id = H5Dopen2(detector_id, "beam_center_x", H5P_DEFAULT);
    H5Dread(beamX_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &header.beam_center[0]);
    H5Dclose(beamX_id);
    beamY_id = H5Dopen2(detector_id, "beam_center_y", H5P_DEFAULT);
    H5Dread(beamY_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &header.beam_center[1]);
    H5Dclose(beamY_id);
    psize_id = H5Dopen2(detector_id, "x_pixel_size", H5P_DEFAULT);
    H5Dread(psize_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &header.psize);
    H5Dclose(psize_id);
    
    // Read Pixel Mask if available
    img.mask = NULL;
    if (H5Lexists(file_id, "/entry/instrument/detector/detectorSpecific/pixel_mask", H5P_DEFAULT) > 0) {
        hid_t mask_id = H5Dopen(file_id, "/entry/instrument/detector/detectorSpecific/pixel_mask", H5P_DEFAULT);
        if (mask_id >= 0) {
            img.mask = (int*)malloc(header.size[0] * header.size[1] * sizeof(int));
            if (H5Dread(mask_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, img.mask) < 0) {
                printf("Warning: Failed to read pixel mask.\n");
                free(img.mask);
                img.mask = NULL;
            }
            H5Dclose(mask_id);
        }
    }
}

Spot* create_spot(int label, int num_pixels) {
    Spot* spot = malloc(sizeof(Spot));
    spot->label = label;
    spot->num_pixels = num_pixels;
    spot->x_coords = malloc(num_pixels * sizeof(int));
    spot->y_coords = malloc(num_pixels * sizeof(int));
    spot->core_distance = UNDEFINED; // Initialize core_distance
    spot->reachability_distance = UNDEFINED; // Initialize reachability_distance
    spot->cluster_id = UNDEFINED; // Initialize cluster_id
    return spot;
}

void filespot(const Spot* spot, int x, int y, FILE* file) {
    fprintf(file, "%4d  %5d     %5d,%5d %8.f     %d\n",
            spot->label, spot->num_pixels, y, x, spot->intensity, spot->cluster_id);
}

void drawcircle(double *image, int spot_radius, int width, int height, int center_x, int center_y, unsigned char circle_color) {
    int thickness = 3;
    for (int x = center_x - spot_radius; x <= center_x + spot_radius; x++) {
        for (int y = center_y - spot_radius; y <= center_y + spot_radius; y++) {
            // Calculate the distance from the center pixel
            int dx = x - center_x;
            int dy = y - center_y;
            int distance = dx * dx + dy * dy;

            bool inside_condition = distance <= spot_radius * spot_radius && distance > (spot_radius - thickness) * (spot_radius - thickness);

            // If the pixel is within the circle radius and on the circle's circumference, set it to the circle color
            if (inside_condition && x >= 0 && x < height && y >= 0 && y < width) {
                image[x * width + y] = circle_color;
            }
        }
    }
}

void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

double calculateresolution(double radius) {
    return (header.wavelength / (sin(atan(radius * header.psize / header.distance) / 2) * 2)); 
}

void applysingulargaussian(double data[], int size, double sigma) {
    double* tempData = malloc(size * sizeof(double));
    memcpy(tempData, data, size * sizeof(double));

    double* kernel = malloc((6 * sigma + 1) * sizeof(double));
    double kernelSum = 0.0;

    int i;
    for (i = 0; i < 6 * sigma + 1; i++) {
        double x = i - 3 * sigma;
        kernel[i] = exp(-0.5 * x * x / (2 * (sigma * sigma)));
        kernelSum += kernel[i];
    }

    for (i = 0; i < 6 * sigma + 1; i++) {
        kernel[i] /= kernelSum;
    }

    for (i = 0; i < size; i++) {
        double blurredValue = 0.0;
        int j;
        for (j = 0; j < 6 * sigma + 1; j++) {
            int index = i + j - 3 * sigma;
            if (index >= 0 && index < size) {
                blurredValue += tempData[index] * kernel[j];
            }
        }
        data[i] = blurredValue;
    }

    free(tempData);
    free(kernel);
}


void addringset(PositionSet** head, double startPosition, double endPosition) {
    PositionSet* newPositionSet = (PositionSet*)malloc(sizeof(PositionSet));
    newPositionSet->startPosition = startPosition;
    newPositionSet->endPosition = endPosition;
    newPositionSet->next = *head;
    *head = newPositionSet;
}

double calculateavgintensity(double* image, int* mask, int width, int height, double thresfactor) {
    long long totalPixels = (long long)width * height;
    
    // Robust mean calculation using sigma clipping to ignore Bragg peaks
    double sum = 0;
    long long count = 0;
    double sum_sq = 0;
    
    // Pass 1: Initial Mean/Std (ignoring masked pixels)
    for (int i = 0; i < totalPixels; i++) {
        if (mask && mask[i] != 0) continue; // Skip bad pixels
        // Also skip very high values directly? No, let stats handle it first
        double val = image[i];
        sum += val;
        sum_sq += val * val;
        count++;
    }
    
    if (count == 0) return 0;
    
    double mean = sum / count;
    double variance = (sum_sq / count) - (mean * mean);
    double std_dev = sqrt(variance > 0 ? variance : 0);
    
    // Pass 2: Sigma Clipping ( exclude > mean + 3*std )
    // This removes the massive spots from the background calculation
    double robust_sum = 0;
    long long robust_count = 0;
    double limit = mean + 3.0 * std_dev;
    
    for (int i = 0; i < totalPixels; i++) {
        if (mask && mask[i] != 0) continue;
        
        if (image[i] <= limit) {
            robust_sum += image[i];
            robust_count++;
        }
    }
    
    if (robust_count == 0) return mean * thresfactor;
    
    return (robust_sum / robust_count) * thresfactor;
}

void calculate_radial_profile(double* data, int width, int height, double* radial_mean, double* radial_max, int max_radius) {
    // efficient radial profile calculation
    // Reset radial arrays
    for (int i = 0; i < max_radius; i++) {
        radial_mean[i] = 0;
        radial_max[i] = 0;
    }
    
    int *counts = (int*)calloc(max_radius, sizeof(int));
    
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int dx = x - header.beam_center[0];
            int dy = y - header.beam_center[1];
            int r = (int)sqrt(dx*dx + dy*dy);
            
            if (r < max_radius) {
                double val = data[y * width + x];
                if (val > 0) { // Should we ignore 0/masked?
                    radial_mean[r] += val;
                    counts[r]++;
                    if (val > radial_max[r]) {
                        radial_max[r] = val;
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < max_radius; i++) {
        if (counts[i] > 0) {
            radial_mean[i] /= counts[i];
        }
    }
    
    free(counts);
}


void intensityanalysis(double data[], double max_data[], double* binary_mask, int width, int height, Spot** spots, int* num_spots_ptr) {
    int num_spots = *num_spots_ptr;
    header.resolution = calculateresolution(header.max_distance);
    int newSize = (img.height + img.width) / 2; 
    
    // Improved Ice Detection Logic
    // Ice rings typically appear at specific resolutions: ~3.67, ~2.25, ~1.92 A (Cubic) or ~3.90, ~2.67, ~2.07 A (Hexagonal)
    // We will look for significant local maxima in the radial profile
    
    // header.ice = 0; // Reset ice flag MOVED to applySPORTA to persist pixel mask flag
    
    // Smoothing the radial profile first (simple moving average)
    int window = 5;
    double *smoothed = (double*)malloc(newSize * sizeof(double));
    for(int i=0; i<newSize; i++) {
        double sum = 0;
        int count = 0;
        for(int w=-window; w<=window; w++) {
            if(i+w >= 0 && i+w < newSize) {
                sum += data[i+w];
                count++;
            }
        }
        smoothed[i] = sum/count;
    }
    
    // Robust Background Estimation using iterative sigma clipping
    double bg_mean = 0;
    double bg_std = 0;
    int n_points = 0;
    
    // Initial pass: collect data in relevant range
    double *subset = (double*)malloc(newSize * sizeof(double));
    for(int i=0; i<newSize; i++) {
        double res = calculateresolution(i);
        if(res > 1.5 && res < 5.0 && data[i] > 0) {
             subset[n_points++] = data[i];
        }
    }
    
    // Iterative clipping
    if (n_points > 5) {
        double current_mean, current_std;
        for (int iter = 0; iter < 3; iter++) {
            double sum = 0, sum_sq = 0;
            int count = 0;
            
            // Calculate stats
            if (iter == 0) {
                for(int i=0; i<n_points; i++) {
                    sum += subset[i];
                    sum_sq += subset[i]*subset[i];
                    count++;
                }
            } else {
                for(int i=0; i<n_points; i++) {
                    if (abs(subset[i] - bg_mean) < 3.0 * bg_std) {
                        sum += subset[i];
                        sum_sq += subset[i]*subset[i];
                        count++;
                    }
                }
            }
            
            if (count < 5) break; 
            
            current_mean = sum / count;
            current_std = sqrt((sum_sq / count) - (current_mean * current_mean));
            
            bg_mean = current_mean;
            bg_std = current_std;
        }
    }
    free(subset);

    // Threshold for general peak detection
    // Lower multiplier (e.g. 2.0 or 2.5) because we stripped outliers (Bragg peaks)
    double threshold = bg_mean + 3.0 * bg_std; 
    if (threshold < bg_mean * 1.2) threshold = bg_mean * 1.2; 
 

    // double start = -1; // Removed unused variable
    double position;

    // Detect peaks using the smoothed profile
    for (int i = 1; i < newSize - 1; i++) {
        position = calculateresolution(i);
        
        // Check if current point is a local maximum and above threshold
        // And within reasonable ice resolution range
        // bool is_peak = (smoothed[i] > smoothed[i-1] && smoothed[i] > smoothed[i+1]); // Removed unused variable
        // bool above_threshold = (smoothed[i] > threshold); // Removed unused variable
        
        // Known ice ring resolutions: 
        // Hexagonal: 3.90, 2.67, 2.07
        // Cubic: 3.67, 2.25, 1.92
        // We define ranges around these values
        bool is_ice_candidate = false;
        if ((position > 3.85 && position < 3.95) || 
            (position > 3.62 && position < 3.72) ||
            (position > 2.62 && position < 2.72) ||
            (position > 2.20 && position < 2.30) ||
            (position > 2.02 && position < 2.12) ||
            (position > 1.87 && position < 1.97)) {
            is_ice_candidate = true;
        }
        
        // Lower threshold for known ice positions
        double local_threshold = is_ice_candidate ? (bg_mean + 2.0 * bg_std) : threshold;
        
        // Also check max profile for "stretchy spots" / arcs
        // If max intensity at this radius is significantly higher than mean, it suggests non-uniformity (spots/arcs)
        // bool strong_max = (max_data[i] > local_threshold * 1.5); // Removed unused variable // Heuristic

    // STRICTER LOGIC: Only flag if it matches known ice resolutions
    // Bragg spots from sample can appear anywhere, we don't want to flag them as ice unless they overlap with ice rings.
    
        if (is_ice_candidate) {
             // For radial profile, only flag if the MEAN is high, indicating a ring or arc.
             // Single spots (high max, low mean) are handled by spot analysis below.
             if (smoothed[i] > local_threshold) {
                 if(header.ice == 0) header.ice = 1; 
                 
                 double width_est = 0.05; 
                 addringset(&positionSets, position - width_est, position + width_est);
                 
                 // Skip neighbors to avoid multiple rings for same peak
                 i += 2; 
             }
        }
    }

    // Detect Ice Rings based on Spot Concentration (for weak rings that appear as spots)
    // Sometimes intensity profile is weak but spots cluster on the ring.
    int ice_counts[6] = {0}; 
    // Ranges centers: 3.90, 3.67, 2.67, 2.25, 2.07, 1.92
    
    for (int k = 0; k < num_spots; k++) {
        Spot* s = spots[k];
        long long sumX = 0, sumY = 0;
        for(int p=0; p<s->num_pixels; p++) {
            sumX += s->x_coords[p];
            sumY += s->y_coords[p];
        }
        double avgX = (double)sumX / s->num_pixels; 
        double avgY = (double)sumY / s->num_pixels; 
        
        double dx = avgY - header.beam_center[0];
        double dy = avgX - header.beam_center[1];
        double r = sqrt(dx*dx + dy*dy);
        double pos = calculateresolution(r);
        
        // Check ranges (approx +/- 0.05 around known values)
        if (pos > 3.85 && pos < 3.95) ice_counts[0]++;
        else if (pos > 3.62 && pos < 3.72) ice_counts[1]++;
        else if (pos > 2.62 && pos < 2.75) ice_counts[2]++; // slightly wider for 2.67
        else if (pos > 2.20 && pos < 2.30) ice_counts[3]++;
        else if (pos > 2.02 && pos < 2.12) ice_counts[4]++;
        else if (pos > 1.87 && pos < 1.97) ice_counts[5]++;
    }
    
    double centers[6] = {3.90, 3.67, 2.67, 2.25, 2.07, 1.92};
    for(int i=0; i<6; i++) {
        if (ice_counts[i] > 10) { // Threshold: >10 spots indicates an ice ring/arc
             if(header.ice == 0) header.ice = 1; 
             // Add ring to set for filtering
             addringset(&positionSets, centers[i] - 0.06, centers[i] + 0.06);
        }
    }
    
    
    // FILTER SPOTS ON ICE RINGS
    if (positionSets != NULL) {
        int kept_spots = 0;
        for (int k = 0; k < num_spots; k++) {
            Spot* s = spots[k];
            // Calculate resolution of spot center
            long long sumX = 0, sumY = 0; // Use long long to prevent overflow
            for(int p=0; p<s->num_pixels; p++) {
                sumX += s->x_coords[p];
                sumY += s->y_coords[p];
            }
            double avgX = (double)sumX / s->num_pixels; 
            double avgY = (double)sumY / s->num_pixels; 
            
            // Distance from beam center
            // double dist = sqrt(pow(avgX - header.beam_center[1], 2) + pow(avgY - header.beam_center[0], 2));  
            // Note: In calculate_radial_profile: dx = x(col) - beam[0], dy=y(row) - beam[1].
            // x_coords/y_coords in spot are row/col?
            // In connectedcomponentanalysis: dfs(i, j...) where i is height(row), j is width(col).
            // dfs stores x_coords[label_count] = i; y_coords[...] = j;
            // So x_coords is ROW (height), y_coords is COL (width).
            // header.beam_center[0] is usually X (column?), [1] is Y (row?). 
            // Check readhdf5: 
            // beamX_id = "beam_center_x" -> header.beam_center[0]
            // beamY_id = "beam_center_y" -> header.beam_center[1]
            // So [0] is col (x), [1] is row (y).
            // Current code in calculate_radial_profile: dx = x - header.beam_center[0]; dy = y - header.beam_center[1];
            // correct.
            // So for Spot: avgX is ROW, avgY is COL.
            // dx = avgY - header.beam_center[0]; 
            // dy = avgX - header.beam_center[1];
            double dx = avgY - header.beam_center[0];
            double dy = avgX - header.beam_center[1];
            
            // Moved logic to end of function
            double r = sqrt(dx*dx + dy*dy);
            
            double resolution = calculateresolution(r);
            
            bool on_ice_ring = false;
            PositionSet* currentSet = positionSets;
            while (currentSet != NULL) {
                // Add a small buffer to the range? The range in positionSet is already +/- 0.05
                if (resolution >= currentSet->startPosition && resolution <= currentSet->endPosition) {
                    on_ice_ring = true;
                    break;
                }
                currentSet = currentSet->next;
            }
            
            if (!on_ice_ring) {
                // Hard Filter for standard ice rings (even if not detected by histogram)
                // Widen tolerance SIGNIFICANTLY to catch "blobs" as reported by user
                // Standard: 3.90, 3.67, 2.67, 2.25, 2.07, 1.92
                // Added: 3.44, 1.72, 1.52 (common ice rings)
                 if ((resolution > 3.80 && resolution < 4.00) || 
                    (resolution > 3.57 && resolution < 3.77) ||
                    (resolution > 3.34 && resolution < 3.54) || // 3.44
                    (resolution > 2.57 && resolution < 2.77) ||
                    (resolution > 2.15 && resolution < 2.35) ||
                    (resolution > 1.97 && resolution < 2.17) ||
                    (resolution > 1.82 && resolution < 2.02) ||
                    (resolution > 1.62 && resolution < 1.82) || // 1.72
                    (resolution > 1.42 && resolution < 1.62)) { // 1.52
                    on_ice_ring = true;
                    if(header.ice == 0) header.ice = 1; 
                }
            }

            if (!on_ice_ring) {
                // Keep this spot
                spots[kept_spots] = s;
                // Debug: Print resolution of kept spots to see what's slipping through
                // printf("Kept Spot: %.3f A (r=%.1f)\n", resolution, r);
                kept_spots++;
            } else {
                // Reject this spot
                // printf("Rejected Spot on Ice Ring: %.3f A\n", resolution);
                
                // Also clear the pixels in the binary mask so it doesn't show up in the output image!
                if (binary_mask != NULL) {
                    for (int p = 0; p < s->num_pixels; p++) {
                        int px = s->x_coords[p];
                        int py = s->y_coords[p];
                        if (px >= 0 && px < height && py >= 0 && py < width) {
                            // px is row, py is col. image is row-major: row * width + col
                            binary_mask[px * width + py] = 0.0;
                        }
                    }
                }

                free(s->x_coords);
                free(s->y_coords);
                free(s);
            }
        }
        *num_spots_ptr = kept_spots;
        num_spots = kept_spots; // Update local var for subsequent loops
        header.num_spots = kept_spots; // Update header too
    }

    // Double-pass: Hard Filter for standard ice rings (even if not detected by histogram)
    // This catches any spots/fragments that survived pixel masking (e.g. if center shifted)
    {
        int kept_spots = 0;
        for (int k = 0; k < num_spots; k++) {
            Spot* s = spots[k];
             double dx = s->y_coords[s->num_pixels/2] - header.beam_center[0];
             double dy = s->x_coords[s->num_pixels/2] - header.beam_center[1];
             double r = sqrt(dx*dx + dy*dy);
             double resolution = calculateresolution(r);
             
             bool on_ice_ring = false;
             if ((resolution > 3.80 && resolution < 4.00) || // 3.90
                 (resolution > 3.57 && resolution < 3.77) || // 3.67
                 (resolution > 3.34 && resolution < 3.54) || // 3.44
                 (resolution > 2.57 && resolution < 2.77) || // 2.67
                 (resolution > 2.15 && resolution < 2.35) || // 2.25
                 (resolution > 1.97 && resolution < 2.17) || // 2.07
                 (resolution > 1.82 && resolution < 2.05) || // 1.95, 1.92, 1.88
                 (resolution > 1.78 && resolution < 2.17) || // Overlap safety
                 (resolution > 1.62 && resolution < 1.82) || // 1.72
                 (resolution > 1.42 && resolution < 1.62)) { // 1.52
                 on_ice_ring = true;
                 if(header.ice == 0) header.ice = 1;
             }
             
             if (!on_ice_ring) {
                 spots[kept_spots] = s;
                 kept_spots++;
             } else {
                 free(s->x_coords);
                 free(s->y_coords);
                 free(s);
             }
        }
        *num_spots_ptr = kept_spots;
        num_spots = kept_spots;
        header.num_spots = kept_spots;
    }
    

    // Debug output status
    // printf("DEBUG: header.ice=%d, positionSets=%p\n", header.ice, (void*)positionSets);


    // Analyze detected spots for "stretchy" ice spots
    for (int k = 0; k < num_spots; k++) {
        Spot* s = spots[k];
        double sumX = 0, sumY = 0;
        int min_x = 10000, max_x = 0;
        int min_y = 10000, max_y = 0;
        
        for(int p=0; p<s->num_pixels; p++) {
            sumX += s->x_coords[p];
            sumY += s->y_coords[p];
            if(s->x_coords[p] < min_x) min_x = s->x_coords[p];
            if(s->x_coords[p] > max_x) max_x = s->x_coords[p];
            if(s->y_coords[p] < min_y) min_y = s->y_coords[p];
            if(s->y_coords[p] > max_y) max_y = s->y_coords[p];
        }
        
        double avgX = sumX / s->num_pixels;
        double avgY = sumY / s->num_pixels;
        
        // Calculate resolution of spot center
        // double dx = avgX - header.beam_center[1]; // Removed unused variable // x is row/height? verify coordinate system. header has beam_center in [0]=X(row?), [1]=Y(col?). 
        // In drawcircle: image[x*width + y]. x is row. header.beam_center[1] is passed as center_x.
        // HDF5 header usually: [0]=X, [1]=Y.
        // main code: dist = sqrt(pow(header.beam_center[0] - j, 2) + pow(header.beam_center[1] - i, 2)); i=row, j=col.
        // So [0] is col (j), [1] is row (i)? 
        // Let's rely on distance calculation from connectedcomponentanalysis:
        // dist = sqrt(pow(header.beam_center[0] - j, 2) + pow(header.beam_center[1] - i, 2));
        // So [0] <-> j (width/col), [1] <-> i (height/row).
        
        // s->x_coords are from `dfs(x, y...) -> pixels[...] = x`. x was i (row).
        // s->y_coords are y (col).
        
        double dist = sqrt(pow(header.beam_center[1] - avgX, 2) + pow(header.beam_center[0] - avgY, 2));
        double position = calculateresolution(dist);
        
        // Check if resolution matches ice
         bool is_ice_res = false;
         if ((position > 3.85 && position < 3.95) || 
            (position > 3.62 && position < 3.72) ||
            (position > 2.62 && position < 2.75) || // widened slightly
            (position > 2.20 && position < 2.30) ||
            (position > 2.02 && position < 2.12) ||
            (position > 1.87 && position < 1.97)) {
            is_ice_res = true;
        }
        
        if (is_ice_res) {
            // Check shape
            double width_box = max_y - min_y + 1;
            double height_box = max_x - min_x + 1;
            double aspect = (width_box > height_box) ? width_box / height_box : height_box / width_box;
            
            // Criteria for "stretchy" or large spot
            // Normal spots are roughly circular (aspect ~1.0-1.5) and small.
            // Stretchy spots: aspect > 2.0?
            // Large spots: num_pixels > 60 (we increased limit to 2000)
            
            // Require minimum size for aspect ratio check to avoid 3x1 pixel noise
            if (s->num_pixels > 60 || (s->num_pixels > 10 && aspect > 1.8)) {
                if(header.ice == 0) header.ice = 1;
                 double width_est = 0.05; 
                 addringset(&positionSets, position - width_est, position + width_est);
            }
        }
    }
    
    free(smoothed);    

    if (cfg.gnu){
        FILE* data_file = fopen("histogram_data.txt", "w");
         if (data_file != NULL) {
            for (int i = 0; i < newSize; i++) {
                double pos = calculateresolution(i);
                if (pos < 10.0) { // Limit plot range
                    fprintf(data_file, "%f %f\n", pos, data[i]);
                }
            }
            fclose(data_file);

            FILE* gnuplot_pipe = popen("gnuplot -persist", "w");
            if (gnuplot_pipe != NULL) {
                fprintf(gnuplot_pipe, "set terminal qt title '%s'\n", cfg.filename);
                fprintf(gnuplot_pipe, "set title 'Radial Profile'\n");
                fprintf(gnuplot_pipe, "set xlabel 'Resolution [Å]'\n");
                fprintf(gnuplot_pipe, "set ylabel 'Intensity'\n");
                fprintf(gnuplot_pipe, "plot 'histogram_data.txt' with lines title 'Radial Profile', %f title 'Threshold'\n", threshold);
                fprintf(gnuplot_pipe, "exit\n");
                pclose(gnuplot_pipe);
            }
            remove("histogram_data.txt");
        }
    }
}


// Helper function to clear pixels in known ice ring resolutions AND beam center
void apply_pixel_masks(double* binary_mask, int width, int height) {
    // Iterate over all pixels
    // 1. Ice Rings: If resolution matches ice ring, set binary_mask to 0
    // 2. Beam Center: If pixel is within 70 pixels of beam center, set to 0
    
    // Beam Center radius squared for fast comparison
    double beam_radius_sq = 70.0 * 70.0;

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            // Calculate distance from beam center
             double dx = j - header.beam_center[0];
             double dy = i - header.beam_center[1];
             double r_sq = dx*dx + dy*dy;
             
             // 1. Beam Center Mask
             if (r_sq < beam_radius_sq) {
                 binary_mask[i*width + j] = 0.0;
                 continue; // Skip ice check if already masked
             }

             double r = sqrt(r_sq);
             double resolution = calculateresolution(r);
             
             // 2. Ice Ring Mask
             // Check against known ice rings with +/- 0.10 A tolerance
             // List: 3.90, 3.67, 3.44, 2.67, 2.25, 2.07, 1.95, 1.92, 1.88, 1.72, 1.52
             if ((resolution > 3.80 && resolution < 4.00) || // 3.90
                 (resolution > 3.57 && resolution < 3.77) || // 3.67
                 (resolution > 3.34 && resolution < 3.54) || // 3.44
                 (resolution > 2.57 && resolution < 2.77) || // 2.67
                 (resolution > 2.15 && resolution < 2.35) || // 2.25
                 (resolution > 1.97 && resolution < 2.17) || // 2.07
                 (resolution > 1.82 && resolution < 2.05) || // 1.95, 1.92, 1.88
                 (resolution > 1.78 && resolution < 2.17) || // Overlap safety
                 (resolution > 1.62 && resolution < 1.82) || // 1.72
                 (resolution > 1.42 && resolution < 1.62)) { // 1.52
                 
                 // Mark pixel as 0
                 binary_mask[i*width + j] = 0.0;
                 if(header.ice == 0) header.ice = 1;
             }
        }
    }
}


void dfs(int x, int y, int current_label, int height, int width, double *image, int* num_pixels, int* pixels, int limit) {
    if (x < 0 || x >= height || y < 0 || y >= width || image[x * width + y] != 1 || *num_pixels > limit) {
        return;
    }

    // Label the current pixel
    image[x * width + y] = current_label;

    // Store the pixel coordinates
    pixels[(*num_pixels) * 2] = x;
    pixels[(*num_pixels) * 2 + 1] = y;
    (*num_pixels)++;

    // Visit the neighbors
    for (int direction = 0; direction < 4; direction++) {
        dfs(x + dx[direction], y + dy[direction], current_label, height, width, image, num_pixels, pixels, limit);
    }
}

Spot** connectedcomponentanalysis(int height, int width, double *image, double *data, double *preradial, double *postradial, int* num_spots) {
    int next_label = 2; // Start at 2 to avoid conflict with 1 (binary mask) and 0 (background)
    int min_pixel = cfg.min_pixel_size;
    int max_pixel = 2000; // Increased to allow for large "stretchy" ice spots
    Spot** spots = NULL;
    int max_spots = 0;
    // int total_blobs = 0; // Removed unused variable
    // int total_blobs = 0; // Removed unused variable
    // int rejected_small = 0; // Removed unused variable
    // int rejected_large = 0; // Removed unused variable
    // int rejected_scd = 0; // Removed unused variable
    // int accepted = 0; // Removed unused variable

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (image[i * width + j] == 1) {
                // total_blobs++; // Removed unused variable
                int distance = (int) sqrt(pow(header.beam_center[0] - j, 2) + pow(header.beam_center[1] - i, 2));
                int current_label = next_label;
                int num_pixels = 0;
                double intensity = 0;
                int* pixels = malloc(height * width * 2 * sizeof(double));
                bool scdcondition = false;
                int encounter = 0;
                dfs(i, j, current_label, height, width, image, &num_pixels, pixels, max_pixel);
                for (int q = -cfg.proximity; q <= cfg.proximity; q++) {
                    for (int w = -cfg.proximity; w <= cfg.proximity; w++) {
                        int checkX = q + i;
                        int checkY = w + j;
                        if (checkX >= 0 && checkX < width && checkY >= 0 && checkY < height) {
                            if (image[checkX * width + checkY] != 0 && image[checkX * width + checkY] != current_label) {
                                encounter++;
                                if (encounter > 0) {
                                    scdcondition = true;
                                    break;
                                }
                            }
                        }
                    }
                    if (scdcondition) {
                        break;
                    }
                }
                
                // if (scdcondition) rejected_scd++;
                // else if (num_pixels < min_pixel) rejected_small++;
                // else if (num_pixels > max_pixel) rejected_large++;
                // else accepted++;
                
                // Allow neighboring spots (disable scdcondition rejection)
                // if (!scdcondition && num_pixels >= min_pixel && num_pixels <= max_pixel) {
                if (num_pixels >= min_pixel && num_pixels <= max_pixel) {
                    if (postradial[distance] == 0){
                        postradial[distance] = data[i * width + j];
                    }
                    else {
                        postradial[distance] += data[i * width + j];
                    }
                    if (*num_spots >= max_spots) {
                        max_spots += 10;
                        spots = realloc(spots, max_spots * sizeof(Spot*));
                    }
                    spots[*num_spots] = create_spot(current_label, num_pixels);
                    spots[*num_spots]->x_coords = malloc(num_pixels * sizeof(double));
                    spots[*num_spots]->y_coords = malloc(num_pixels * sizeof(double));
                    for (int k = 0; k < num_pixels; k++) {
                        spots[*num_spots]->x_coords[k] = pixels[k * 2];
                        spots[*num_spots]->y_coords[k] = pixels[k * 2 + 1];
                        intensity += data[pixels[k * 2] * width + pixels[k * 2 + 1]];
                    }
                    spots[*num_spots]->intensity = intensity;
                    (*num_spots)++; 
                }
                
                free(pixels);
                next_label++;
            }
        }
    }
    // printf("\nDEBUG Spot Detection Stats:\n");
    // printf("Total blobs found: %d\n", total_blobs);
    // printf("Rejected (Small < %d): %d\n", min_pixel, rejected_small);
    // printf("Rejected (Large > %d): %d\n", max_pixel, rejected_large);
    // printf("Rejected (Neighbor Conflict SCD): %d\n", rejected_scd);
    // printf("Accepted: %d\n\n", accepted);
    return spots;
}

int comparedistances(const void* a, const void* b) {
    double distanceA = *(const double*)a;
    double distanceB = *(const double*)b;

    if (distanceA < distanceB) {
        return -1;
    } else if (distanceA > distanceB) {
        return 1;
    } else {
        return 0;
    }
}

double euclideandistance(Spot* a, Spot* b) {
    double dx = a->x_coords[0] - b->x_coords[0];
    double dy = a->y_coords[0] - b->y_coords[0];
    return sqrt(dx * dx + dy * dy);
}

void updatecoredistance(Spot* spot, Spot** neighbors, int num_neighbors, double epsilon, unsigned int minPts) {
    if (num_neighbors >= minPts) {
        double* distances = (double*)malloc(num_neighbors * sizeof(double));
        for (int i = 0; i < num_neighbors; ++i) {
            distances[i] = euclideandistance(spot, neighbors[i]);
        }
        qsort(distances, num_neighbors, sizeof(double), comparedistances);

        spot->core_distance = distances[minPts - 1] + epsilon;

        free(distances);
    } else {
        spot->core_distance = UNDEFINED;
    }
}

void optics(Spot** spots, int num_spots, double epsilon, unsigned int minPts) {
    for (int i = 0; i < num_spots; ++i) {
        if (spots[i]->cluster_id != UNCLASSIFIED) {
            continue;
        }

        Spot** neighbors = (Spot**)malloc(num_spots * sizeof(Spot*));
        int num_neighbors = 0;

        for (int j = 0; j < num_spots; ++j) {
            if (i == j) {
                continue;
            }

            double distance = euclideandistance(spots[i], spots[j]);
            if (distance <= epsilon) {
                neighbors[num_neighbors] = spots[j];
                ++num_neighbors;
            }
        }

        spots[i]->cluster_id = NOISE;
        spots[i]->reachability_distance = UNDEFINED;

        if (num_neighbors >= minPts) {
            updatecoredistance(spots[i], neighbors, num_neighbors, epsilon, minPts);

            for (int j = 0; j < num_neighbors; ++j) {
                Spot* neighbor = neighbors[j];
                if (neighbor->cluster_id == UNCLASSIFIED) {
                    double reachability_distance = euclideandistance(spots[i], neighbor);
                    if (reachability_distance <= epsilon) {
                        neighbor->reachability_distance = fmax(reachability_distance, spots[i]->core_distance);
                    }
                }
            }

            spots[i]->cluster_id = CORE_POINT;
        }

        free(neighbors);
    }
}

// Comparison function for qsort (Descending Angstroms: High -> Low)
// So index 0 = worst resolution (large A), index 95% = near best (small A)
int compare_doubles(const void* a, const void* b) {
    double arg1 = *(const double*)a;
    double arg2 = *(const double*)b;
    if (arg1 < arg2) return 1; // Descending
    if (arg1 > arg2) return -1;
    return 0;
}

void plot3DImage(Spot** spots, int width, int height) {
    FILE* data_file = fopen("image_data.txt", "w");
    if (data_file == NULL) {
        fprintf(stderr, "Error opening data file\n");
        return;
    }

    for (int i = 0; i < header.num_spots; i++) {
        int index = spots[i]->num_pixels;
        for (int k = 0; k < index; k++) {
            int cluster_id = spots[i]->cluster_id;
            fprintf(data_file, "%d %d %d\n", spots[i]->x_coords[k], spots[i]->y_coords[k], cluster_id);
        }
    }
    fclose(data_file);

    // Plot the 3D data using gnuplot
    FILE* gnuplot_pipe = popen("gnuplot -persist", "w");
    if (gnuplot_pipe == NULL) {
        fprintf(stderr, "Error opening Gnuplot pipe\n");
        return;
    }

    fprintf(gnuplot_pipe, "set terminal qt title '3D Cluster Plot'\n");
    // fprintf(gnuplot_pipe, "set xlabel 'X' font ',12'\n");
    // fprintf(gnuplot_pipe, "set ylabel 'Y' font ',12'\n");
    // fprintf(gnuplot_pipe, "set zlabel 'Cluster ID' font ',12'\n");
    fprintf(gnuplot_pipe, "unset key\n");
    fprintf(gnuplot_pipe, "set xtics font ',9'\n");
    fprintf(gnuplot_pipe, "set ytics font ',9'\n");
    fprintf(gnuplot_pipe, "unset ztics font ',11'\n");

    // Specify custom binary palette for cluster IDs
    fprintf(gnuplot_pipe, "set palette defined (0 'dark-red', 1 'dark-blue')\n");

    // Remove the color scale legend (color box)
    fprintf(gnuplot_pipe, "unset colorbox\n");

    // Plot the data points in 3D space with colors representing the cluster_id
    fprintf(gnuplot_pipe, "plot 'image_data.txt' using 1:2:3 with points palette pointtype 7 pointsize 0.4\n");
    fprintf(gnuplot_pipe, "exit\n");

    pclose(gnuplot_pipe);
}



void applySPORTA(double* image, double* data, int height, int width){
    int num_spots = 0;
    header.ice = 0; // Initialize ice flag here
    double radial[img.height + img.width], postradial[img.height + img.width];
    
    // We remove the old binary mask intensity calculation that resulted in 0.00
    // header.postintensity = calculateavgintensity(img.data, width, height, 1);
    
    // Clear ice rings and beam center from binary mask BEFORE spot detection
    // This ensures "streaky" spots or blobs on ice rings are never detected in the first place
    apply_pixel_masks(image, width, height);
    
    Spot** spots = connectedcomponentanalysis(height, width, image, data, radial, postradial, &num_spots);
    header.num_spots = num_spots;

    double epsilon = 90;
    int minPts = 2;
    optics(spots, num_spots, epsilon, minPts);

    unsigned char circle_color = 255;
    FILE* file = fopen("output.txt", "w");
    fprintf(file, "Label  Pixels    Coordinate   Intensity   ClusterID\n");
    // double pX = 0;
    // double tightness = 100;
    // double distance = 0;
    // double max_distance = 0;

    // Spot* furthest_spot = NULL;
    double max_distance_from_center = 0;

    double totalintensity = 0; // Use double for intensity

    // Array to store resolutions for 95% calculation
    double* spot_resolutions = (double*)malloc(num_spots * sizeof(double));
    int res_count = 0;

    for (int i = 0; i < num_spots; i++) {
        int totalX = 0;
        int totalY = 0;
        int num = spots[i]->num_pixels;
        for (int k = 0; k < num; k++) {
            totalX += spots[i]->x_coords[k];
            totalY += spots[i]->y_coords[k];
        }
        
        // Sum intensity for this spot
        totalintensity += spots[i]->intensity;
        
        int aX = (int)(totalX / num);
        int aY = (int)(totalY / num);
        if (spots[i]->cluster_id == 0) {
            // double spot_distance_from_center = sqrt(pow(header.beam_center[1] - aX, 2) + pow(header.beam_center[0] - aY, 2));
            // if (spot_distance_from_center > max_distance_from_center) {
            //    max_distance_from_center = spot_distance_from_center;
                // furthest_spot = spots[i];
            // }
        }
        
        // Calculate distance for ALL spots regardless of cluster_id
        double spot_distance_from_center = sqrt(pow(header.beam_center[1] - aX, 2) + pow(header.beam_center[0] - aY, 2));
        if (spot_distance_from_center > max_distance_from_center) {
            max_distance_from_center = spot_distance_from_center;
        }

        // Store resolution for percentile calculation
        // We can reuse the `radial` array momentarily or alloc new one? 
        // Alloc new one is safer.
        // Wait, we need to collect ALL resolutions to sort.
        // Let's do it after the loop if we store them.. but we don't store them in spots struct?
        // Spots struct has x_coords/y_coords arrays.
        // We can calculate resolution here and store in a dynamic array.


        if (cfg.output){
            filespot(spots[i], aX, aY, file);
            drawcircle(image, 10, width, height, aX, aY, circle_color);
        }
        
        // Store resolution
        // Calculate resolution from distance
        // r = distance
        // resolution = lambda / (2 * sin(0.5 * atan(r * pixel_size / distance)))
        // Use existing calculateresolution function? Yes.
        double r = spot_distance_from_center;
        if (spot_resolutions) {
            spot_resolutions[res_count++] = calculateresolution(r);
        }
    }
    
    // Calculate 95th Percentile Resolution
    // Sort resolutions descending (High Å -> Low Å). Index 95% = near-best resolution.
    header.resolution_95 = 0.0;
    if (spot_resolutions && res_count > 0) {
        qsort(spot_resolutions, res_count, sizeof(double), compare_doubles);
        int idx95 = (int)(0.95 * res_count);
        if (idx95 >= res_count) idx95 = res_count - 1;
        header.resolution_95 = spot_resolutions[idx95];
        free(spot_resolutions);
    }
    
    // Update header postintensity to be Avg Spot Intensity
    if(num_spots > 0) {
        header.postintensity = totalintensity / num_spots;
    } else {
        header.postintensity = 0.0;
    }
    
    // Set header.intensity to Total Integrated Intensity
    header.intensity = totalintensity;
    
    header.snr = log10(header.num_spots / img.preintensity);
    // Adjusted max_distance calculation without 0.95 factor
    header.max_distance = max_distance_from_center; 
    if (header.max_distance <= 0) header.max_distance = 1.0; // Avoid potential issues with 0 radius

    if (cfg.output){
        drawcircle(image, header.max_distance, width, height, header.beam_center[1], header.beam_center[0], circle_color);
    }
    
    // Calculate Radial Profile
    // int max_radius = (img.height > img.width) ? img.height : img.width; // Removed unused variable 
    // Usually diagonals are larger, but let's be safe with buffer size equal to max dim or diagonal
    // The previous code seemed to use (h+w)/2 which is odd, let's stick to diagonal for safety
    int diagonal = (int)sqrt(width*width + height*height) + 1;
    // Re-allocate localized radial arrays to ensure size correctness
    // Note: The original code used variable length array on stack (VLA) which can stack overflow for large images
    // double radial[img.height + img.width] is risky.
    // Let's use heap allocation.
    
    double *safe_radial = (double*)malloc(diagonal * sizeof(double));
    double *safe_radial_max = (double*)malloc(diagonal * sizeof(double));
    
    if (safe_radial && safe_radial_max) {
        calculate_radial_profile(data, width, height, safe_radial, safe_radial_max, diagonal);
        // Pass 'image' (which is img.data, the binary mask) to be updated
        intensityanalysis(safe_radial, safe_radial_max, image, width, height, spots, &num_spots);
        free(safe_radial);
        free(safe_radial_max);
    } else {
        // Fallback or error
        if(safe_radial) free(safe_radial);
        if(safe_radial_max) free(safe_radial_max);
        printf("Error allocating memory for radial profile.\n");
    }
    // plot3DImage(spots, width, height);
    fclose(file); 
    free(spots);
}


void applythreshold(double* image, double* data, int* mask, int width, int height) {
    // Note: Histogram calculation for Otsu might need update to ignore mask, 
    // but here we use 'threshold' from header which is set by this function?
    // Wait, applythreshold calculates the threshold using Otsu logic on 'image' (which is blurred)
    // The previous logic was:
    // ... calculate Otsu ...
    // then override threshold with header.threshold = threshold * thresfactor
    // But 'calculateavgintensity' sets 'img.preintensity' which is used for SNR, 
    // it does NOT set 'header.threshold'.
    // Actually, 'applythreshold' calculates the threshold.
    // AND 'calculateavgintensity' is called BEFORE 'applythreshold'.
    
    // The user wants 'weak spots'.
    // Otsu on diffraction images is tricky because background > foreground pixels.
    // The previous 'calculateavgintensity' returned a value.
    // Let's use THAT value as the threshold base instead of Otsu?
    // Or improve Otsu?
    // The code currently does:
    
    // 1. img.preintensity = calculateavgintensity(...) -> robust mean
    // 2. applygaussianblur
    // 3. applythreshold -> calculates Otsu on BLURRED image, then sets header.threshold.
    
    // User complaint: "image background is set too high".
    // This implies the threshold calculated HERE is too high.
    // Otsu tries to separate background from foreground (spots).
    // If spots are sparse, Otsu might fail or set threshold too high if bimodality is weak?
    // actually, for sparse spots, Otsu often sets threshold too LOW (just above background noise).
    // If the user says it's too high, maybe the "background" peak is wide?
    
    // PROPOSAL:
    // Use the robust mean from 'calculateavgintensity' as the threshold base!
    // Simply: threshold = img.preintensity (which already includes thresfactor).
    // This gives direct control via -t flag and uses the robust background.
    
    double threshold = img.preintensity; // Already has thresfactor applied
    header.threshold = threshold;

    for (int i = 0; i < height * width; ++i) {
        // Apply Mask - DISABLED to see if it removes valid saturated spots
        // if (mask && mask[i] != 0) {
        //    image[i] = 0.0;
        //    continue;
        // }
        
        int pixelValue = (int)image[i];

        // If the pixel intensity is greater than the threshold
        if (pixelValue > threshold) {
            image[i] = 1.0;
        }
        // Otherwise, set it to 0
        else {
            image[i] = 0.0;
        }
    }
}


void applygaussianblur(double* image, int width, int height) {
    double sigma;
    if (cfg.gaussian > 0) sigma = cfg.gaussian;
    else{
        return;
    }
    sigma = 0.1;
    // Create the Gaussian kernel
    int kernel_size = ceil(6 * sigma) + 1;
    int radius = kernel_size /2;
    float* kernel = (float*)malloc(kernel_size * sizeof(float));
    float kernel_sum = 0.0;

    for (int i = 0; i < kernel_size; i++) {
        int x = i - radius;
        kernel[i] = exp(-(x * x) / (3 * sigma * sigma));
        kernel_sum += kernel[i];
    }

    // Normalize the kernel
    for (int i = 0; i < kernel_size; i++) {
        kernel[i] /= kernel_sum;
    }

    // Create a temporary buffer to store the blurred image
    float* blurred_image = (float*)malloc(width * height * sizeof(float));

    // Convolve each pixel with the Gaussian kernel
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            float sum = 0.0;
            int count = 0;

            // Apply the kernel to the current pixel
            for (int i = 0; i < kernel_size; i++) {
                int image_x = x + i - radius;
                if (image_x >= 0 && image_x < width) {
                    sum += image[y * width + image_x] * kernel[i];
                    count++;
                }
            }

            // Store the blurred pixel in the temporary buffer
            blurred_image[y * width + x] = sum / count;
        }
    }

    // Copy the blurred image back to the original image buffer
    for (int i = 0; i < width * height; i++) {
        image[i] = blurred_image[i];
    }

    // Free the memory
    free(kernel);
    free(blurred_image);
}
char* generatemaster(const char* filename) {
    // Check if the input file exists
    if (access(filename, F_OK) == -1) {
        printf("File %s does not exist.\n", filename);
        return NULL;
    }
    
    // If the file is already a master file, return it
    if (strstr(filename, "_master.h5") != NULL) {
        return strdup(filename);
    }

    // Copy the filename to a new variable that we can modify
    static char master[256];
    strncpy(master, filename, sizeof(master));
    master[sizeof(master) - 1] = '\0';  // Ensure null termination

    // Find the last occurrence of "_data_" in the filename
    char *substring = strstr(master, "_data_");
    if (substring == NULL) {
        // Fallback: assume it might be a master file named differently or user error
        // But for now, just print error as before if it doesn't match expected pattern
        printf("Invalid filename format. Expected '_data_' or '_master.h5'\n");
        return NULL;
    }

    // Replace "_data_" and everything after it with "_master.h5"
    memcpy(substring, "_master.h5", strlen("_master.h5"));
    // Make sure the new filename string is properly terminated
    substring[strlen("_master.h5")] = '\0';

    // Check if the new file exists
    if (access(master, F_OK) == -1) {
        printf("File %s does not exist.\n", master);
        return NULL;
    }

    return master;
}

double* deepcopyimage(const double* source, int width, int height) {
    int dataSize = width * height;
    double* target = (double*)malloc(dataSize * sizeof(double));

    for (int i = 0; i < dataSize; i++) {
        target[i] = source[i];
    }

    return target;
}

void setdefault(){
    header.size[0] = 3110;
    header.size[1] = 3269;
    header.distance = 0.2;
    header.wavelength = 0.953742;
    header.psize = 0.000075;
    header.beam_center[0] =1544;
    header.beam_center[1] = 1577;
}

// ends_with_data, contains_proc, contains_native - commented out as unused
// int ends_with_data(const char *path, char *folder) { return 0; }
// int contains_proc(const char *name) { return strstr(name, "proc") != NULL; }
// int contains_native(const char *name) { return strstr(name, "native") != NULL; }

int screener(const char *filename){
    cfg.filename = (char*)filename; // Ensure global is in sync if needed elsewhere
    cfg.master = generatemaster(filename);
    if (cfg.master != NULL){
        readhdf5(cfg.master);
    }
    else{
        setdefault();
    }
    
    char dataset_path[256];
    int read_index = cfg.index; // Default to using configured index directly (if not master/chunked)

    // Robust dataset determination logic
    if (strstr(filename, "_master.h5") != NULL) {
        hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id >= 0) {
            // Priority 1: Check for consolidated 'data' dataset
            if (H5Lexists(file_id, "/entry/data/data", H5P_DEFAULT) > 0) {
                strcpy(dataset_path, "/entry/data/data");
                read_index = cfg.index - 1; // 1-based to 0-based
                if (read_index < 0) read_index = 0;
            } 
            // Priority 2: Check for chunked 'data_000001'
            else if (H5Lexists(file_id, "/entry/data/data_000001", H5P_DEFAULT) > 0) {
                 hid_t dset = H5Dopen(file_id, "/entry/data/data_000001", H5P_DEFAULT);
                 if (dset >= 0) {
                     hid_t space = H5Dget_space(dset);
                     hsize_t dims[3];
                     H5Sget_simple_extent_dims(space, dims, NULL);
                     H5Sclose(space);
                     H5Dclose(dset);
                     
                     int chunk_size = (int)dims[0]; // typically 100 or 1000
                     if (chunk_size > 0) {
                         int eff_index = cfg.index - 1; // 1-based to 0-based
                         if (eff_index < 0) eff_index = 0;
                         
                         int chunk_id = (eff_index / chunk_size) + 1;
                         int local_idx = eff_index % chunk_size;
                         
                         snprintf(dataset_path, sizeof(dataset_path), "/entry/data/data_%06d", chunk_id);
                         read_index = local_idx;
                     } else {
                         // Fallback if dim is 0?
                         snprintf(dataset_path, sizeof(dataset_path), "/entry/data/data_%06d", cfg.index);
                     }
                 } else {
                     snprintf(dataset_path, sizeof(dataset_path), "/entry/data/data_%06d", cfg.index);
                 }
            } else {
                // Fallback: simple numeric mapping (old behavior)
                snprintf(dataset_path, sizeof(dataset_path), "/entry/data/data_%06d", cfg.index);
            }
            H5Fclose(file_id);
        } else {
             // File open failed, but let decompresshdf5 handle error reporting
             snprintf(dataset_path, sizeof(dataset_path), "/entry/data/data_%06d", cfg.index);
        }
    } else {
        strncpy(dataset_path, cfg.data, sizeof(dataset_path));
        // If not master, we assume single file, maybe 1-based index needs offset?
        // Old code used cfg.index directly. Let's keep it but ideally it should be -1 if 1-based.
        // For compatibility with old behavior on single files, we keep cfg.index unless we are sure.
    }
    
    double *data = decompresshdf5(filename, dataset_path, read_index);
    if (data == NULL) {
        fprintf(stderr, "Failed to decompress/read HDF5 data from %s (dataset: %s)\n", filename, dataset_path);
        return 1;
    }
    
    img.height = (int)header.size[0];
    img.width = (int)header.size[1];
    img.data = deepcopyimage(data, img.width, img.height);
    img.preintensity = calculateavgintensity(data, img.mask, img.width, img.height, cfg.thresfactor);
    applygaussianblur(img.data, img.width, img.height);
    applythreshold(img.data, data, img.mask, img.width, img.height);
    applySPORTA(img.data, data, img.width, img.height);
    if (cfg.output){
        savepgm("pre-output", data, img.height, img.width);
        savepgm("post-output", img.data, img.height, img.width);
    }
    mainprint();
    mainprint();
    if(img.mask) free(img.mask);
    free(data);
    free(img.data);

    return 0;
}

int matches_pattern(const char *filename) {
    regex_t regex;
    int ret;
    // Compile the regular expression
    if (regcomp(&regex, ".*_000001\\.h5$", REG_EXTENDED) != 0) {
        printf("Could not compile regex\n");
        return 0;
    }
    // Execute the regular expression
    ret = regexec(&regex, filename, 0, NULL, 0);
    // Free up the regular expression
    regfree(&regex);
    // Check the result
    if (ret == 0) {
        return 1;
    } else {
        return 0;
    }
}

// (getrecentfile, free_memory, getdirs functions removed - unused)


int argparse(int argc, char **argv) {
    optind = 0; // Reset getopt state for re-entrancy (0 is safer on some systems to force re-scan)
    cfg.master = NULL;
    cfg.output = false;
    cfg.index = 1;
    cfg.data = "/entry/data/data";
    cfg.thresfactor = 6;
    cfg.proximity = 2;
    cfg.percent = 0.95;
    cfg.traversing = 0;  
    cfg.searchdir = "data";
    cfg.gnu = false;
    cfg.min_pixel_size = 5; // Default to 5 as requested
    header.ice = 0;
    int c;

    while ((c = getopt(argc, argv, "hof:i:d:t:s:rp:g:m:")) != -1) {
        switch (c) {
            case 'h':
                printf("Usage: ./imager -f <filepath> -o\n");
                printf("Options:\n");
                printf("  -h          Display this help message\n");
                printf("  -f <file>   Specify the input file\n");
                printf("  -o          Enable output\n");
                printf("  -i <index>  Set the index to a specified value\n");
                printf("  -d <data>   Set the data to a specified value\n");
                printf("  -r <dir>    Enable traversing on path <dir>\n");
                printf("  -t <value>  Set the threshold scale factor to a specified value\n");
                printf("  -s <dir>    Set the search directory to a specified value\n");
                printf("  -e <value>  Set the ring exclusion proximity radius\n");
                printf("  -p          Enable histogram through GNUPLOT\n");
                printf("  -g <sigma>  Set gaussian filtering sigma value\n");
                printf("  -m <value>  Set minimum pixel size for spots (default: 5)\n");
                printf("  -n          Display the notes, tips, and comments");
                return(0);
            case 'd':
                cfg.data = optarg;
                break;
            case 'o':
                cfg.output = true;
                break;
            case 'i':
                cfg.index = atoi(optarg);
                break;
            case 'r':
                cfg.traversing =1;
                break;   
            case 't':
                cfg.thresfactor = atoi(optarg);
                break;
            case 's':
                cfg.searchdir =optarg;
                break;
            case 'e':
                cfg.percent =atof(optarg);
                break;
            case 'p':
                cfg.gnu = true;
                break;
            case 'g':
                cfg.gaussian = atof(optarg);
                break;
            case 'f':
                cfg.filename = optarg;
                break;
            case 'm':
                cfg.min_pixel_size = atoi(optarg);
                break;  
            case 'n':
                noteprint();
                return 0;
            case '?':
                if (optopt == 'f' || optopt == 'i' || optopt == 't' || optopt == 's' || optopt == 'd') {
                    fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                } else if (isprint(optopt)) {
                    fprintf(stderr, "Unknown option `-%c'.\n", optopt);
                } else {
                    printf("Help: use -f for file input, -o for output\n");

                }
                return 1;
            default:
                abort();
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    argparse(argc, argv);
    return 0;
}

Hdf5Header getheaderdata() {
    return header;
}

