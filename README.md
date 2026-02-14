# SPORTA

**Spot Detection & Screening of X-Ray Diffraction Images**

SPORTA is a tool for analyzing X-ray diffraction images (HDF5 format) to detect spots, calculate resolution, and screen for ice rings. It combines a high-performance C extension for image processing with a Python interface for ease of use and data handling.

## Features

- fast HDF5 image reading (supports Bitshuffle/LZ4 compression)
- spot detection and connected component analysis
- resolution calculation based on experimental parameters
- ice ring detection
- filtering and thresholding
- batch processing of multiple files

## Installation

SPORTA is designed to be easily installable on Linux environments.

### Prerequisites

- Python 3.x
- HDF5 libraries (installed via system package manager or custom build)
- GCC compiler

### Quick Start

1.  Clone the repository:
    ```bash
    git clone https://github.com/denis/sporta.git
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
sporta /path/to/data/image_001.h5

# Process multiple files using wildcards
sporta /path/to/data/image_*.h5
```

If no arguments are provided, SPORTA defaults to scanning the current directory for a `data` folder and processing files within it that match internal patterns.

### Options

-   `-i <index>`: Set the internal frame index to process within the HDF5 file (default: 1).
    ```bash
    sporta -i 5 image.h5
    ```
-   `-t <value>`: Set the threshold scale factor (default: 6).
-   `-e <value>`: Set the ring exclusion proximity radius (default: 2).
-   `-g <sigma>`: Set the Gaussian filtering sigma value (default: various).
-   `-o`: Enable output of intermediate PGM images (`pre-output.pgm`, `post-output.pgm`).

## Output

The tool prints analysis results to the console and saves a summary to `data.tsv`.

## License

MIT License
