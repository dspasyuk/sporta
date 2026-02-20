# SPORTA

**Spot Detection & Screening of X-Ray Diffraction Images**

SPORTA is a tool for analyzing X-ray diffraction images (HDF5 format) to detect spots, calculate resolution, and screen for ice rings. It combines a high-performance C extension for image processing with a Python interface for ease of use and data handling.

## Features

- Fast HDF5 image reading (supports Bitshuffle/LZ4 compression)
- Spot detection and connected component analysis
- Resolution calculation based on experimental parameters (best spot + 95th percentile)
- Ice ring detection (hexagonal and cubic)
- Gaussian filtering and thresholding
- Minimum pixel size filtering for spots
- Batch processing of multiple files and frame ranges

## Installation

SPORTA is designed to be easily installable on Linux environments.

### Prerequisites

- Python 3.x
- HDF5 libraries (installed via system package manager or custom build)
- GCC compiler

### Quick Start

1.  Clone the repository:
    ```bash
    git clone https://github.com/dspasyuk/sporta.git
    cd sporta
    ```

2.  Run the installation script:
    ```bash
    ./install.sh
    ```

    This script will:
    - Check for system dependencies (and warn if missing).
    - Create a Python virtual environment (`venv`).
    - Install Python dependencies (`pandas`, `numpy`).
    - Compile the C extension and install the `sporta` package.

### Custom HDF5 Path

If you have HDF5 in a custom location (e.g., if the install script warns about missing headers), set the `HDF5_HOME` environment variable before running `install.sh`:

```bash
export HDF5_HOME=/path/to/hdf5-installation
./install.sh
```

## Usage

Activate the virtual environment:

```bash
source venv/bin/activate
```

Then run `sporta`.

### Processing Files

You can process specific HDF5 files by passing them as arguments:

```bash
# Process a single file
sporta /path/to/data/image_master.h5

# Process multiple files using wildcards
sporta /path/to/data/image_*.h5
```

If no arguments are provided, SPORTA defaults to scanning the current directory for a `data` folder and processing files within it that match internal patterns.

### Options

-   `-v`: Verbose mode. Enables saving of analysis results to `data.tsv`.
-   `-i <index>`: Set the internal frame index to process within the HDF5 file (default: 1).
    ```bash
    sporta -i 5 image_master.h5
    ```
-   `-r <start>-<end>`: Process a range of frames (inclusive) within the HDF5 file.
    ```bash
    sporta -r 1-10 image_master.h5
    ```
-   `-t <value>`: Set the threshold scale factor (default: 6).
-   `-e <value>`: Set the percentile value for resolution calculation (default: 0.95).
-   `-g <sigma>`: Set the Gaussian filtering sigma value (optional).
-   `-m <size>`: Set the minimum pixel size for spots (default: 5).
-   `-o`: Enable output of intermediate PGM images (`pre-output.pgm`, `post-output.pgm`).
-   `-h`: Display help message.

## Python API

You can also use SPORTA as a Python library to process files programmatically.

```python
import sporta

# Process a single file
result = sporta.process_file(
    "/path/to/data.h5",
    index=1,
    threshold=6.0,
    min_pixel=5
)
print(result)

# Result is a dictionary containing resolution, spot count, etc.
# {
#     'File': '...',
#     'Resolution': 3.88,
#     '# Spots': 27,
#     'Ice': 0,
#     'SNR': 1.45,
#     ...
# }
```

This allows for easy integration into batch processing pipelines or data analysis scripts. See `example.py` for more usage patterns.


## Output

The tool prints analysis results to the console. If `-v` is used, it also saves a summary to `data.tsv`.

### Metrics Definitions

-   **Best Spot Resolution**: Resolution (in Å) of the farthest detected spot from the beam center.
-   **Resolution (95%)**: Resolution at which 95% of detected spots lie, providing a more representative measure.
-   **Detector Edge Resolution**: Theoretical maximum resolution at the detector edge.
-   **Total Intensity**: Sum of intensities of all pixels within identified spots.
-   **Avg Spot Intensity**: Average intensity per spot (Total Intensity / Spot Count).
-   **Avg Raw Intensity**: Average pixel intensity of the entire raw image (sigma-clipped to exclude Bragg peaks).
-   **Signal Quality (LogSNR)**: Logarithmic signal-to-noise ratio based on spot count vs raw intensity.

## License

MIT License

## Acknowledgments

Special thanks to Aaryan Patel for significant contributions to the development of this code.
