import ctypes
import sys
import os
import pandas as pd

import glob

# Try to find the shared library
lib_path = None

# 1. Look in current directory (local build/dev)
local_lib = glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), "sporta_lib*.so"))
if local_lib:
    lib_path = local_lib[0]

# 2. Look in site-packages (installed via pip)
if lib_path is None:
    try:
        import sysconfig
        # Construct expected name or pattern
        # This is tricky because extension naming varies by platform/python version
        # But we can look relative to where this file is installed
        # If sporta.py is in site-packages, the .so should be next to it
        installed_lib = glob.glob(os.path.join(os.path.dirname(__file__), "sporta_lib*.so"))
        if installed_lib:
            lib_path = installed_lib[0]
            
        # 3. Fallback: manual build name
        if lib_path is None:
             manual_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sporta.so")
             if os.path.exists(manual_path):
                 lib_path = manual_path
                 
    except Exception:
        pass

if lib_path is None:
    # 4. Global fallback
    # Maybe we are running a script that imports sporta, and sporta is installed in site-packages
    # We can try to import the extension module directly to get its file?
    try:
        import importlib.util
        spec = importlib.util.find_spec("sporta_lib")
        if spec and spec.origin:
            lib_path = spec.origin
    except ImportError:
        pass

if lib_path is None:
     raise FileNotFoundError("Could not find sporta_lib shared library. Please install the package with 'pip install .'")

lib = ctypes.CDLL(lib_path)

class Hdf5Header(ctypes.Structure):
    _fields_ = [
        ("distance", ctypes.c_double),
        ("wavelength", ctypes.c_double),
        ("beam_center", ctypes.c_double * 2),
        ("size", ctypes.c_double * 2),
        ("psize", ctypes.c_double),
        ("threshold", ctypes.c_double),
        ("radius_90_percent", ctypes.c_int),
        ("resolution", ctypes.c_double),
        ("num_spots", ctypes.c_int),
        ("intensity", ctypes.c_double),
        ("ice", ctypes.c_int),
        ("max_distance", ctypes.c_double),
        ("snr", ctypes.c_double),
        ("postintensity", ctypes.c_double),
        ("resolution_95", ctypes.c_double),
    ]


lib.getheaderdata.restype = Hdf5Header 

def getdirs(currentdir, folder):
    """Pure-Python directory finder (C version was removed as unused)."""
    import os
    results = []
    for root, dirs, files in os.walk(currentdir):
        if os.path.basename(root) == folder and 'proc' not in root and 'native' not in root:
            results.append(root)
    return results

lib.screener.restype = ctypes.c_int
lib.screener.argtypes = [ctypes.c_char_p]

def process_file(filename, index=1, threshold=6, proximity=2, gaussian=None, min_pixel=5, output=False):
    """
    Process a single HDF5 file or master file.
    
    Args:
        filename (str): Path to the .h5 file.
        index (int): Frame index to process (1-based, default: 1).
        threshold (float): Threshold scale factor (default: 6).
        proximity (float): Ring exclusion proximity radius (default: 2).
        gaussian (float): Gaussian sigma (optional).
        min_pixel (int): Minimum pixel size for spots (default: 5).
        output (bool): Whether to save intermediate PGM images (default: False).
        
    Returns:
        dict: A dictionary containing the analysis results (Resolution, Spot Count, etc.)
    """
    abs_path = os.path.abspath(filename)
    if not os.path.exists(abs_path):
        raise FileNotFoundError(f"File not found: {abs_path}")
        
    # Construct arguments for argparse
    # argparse expects [program_name, arg1, val1, ...]
    args = ["sporta"]
    args.extend(["-t", str(threshold)])
    args.extend(["-e", str(proximity)])
    args.extend(["-i", str(index)])
    if gaussian is not None:
        args.extend(["-g", str(gaussian)])
    args.extend(["-m", str(min_pixel)])
    if output:
        args.append("-o")
        
    # Convert to C-compatible argv
    argc = len(args)
    c_args = (ctypes.c_char_p * (argc + 1))()
    c_args[:-1] = [arg.encode('utf-8') for arg in args]
    c_args[-1] = None
    
    # Update global config in C
    lib.argparse(argc, c_args)
    
    # Run analysis
    res = lib.screener(abs_path.encode())
    if res != 0:
        raise RuntimeError(f"Analysis failed for {abs_path}")
        
    # Retrieve results
    header = lib.getheaderdata()
    return {
        'File': abs_path,
        'Resolution': header.resolution,
        'Distance': header.distance,
        'Wavelength': header.wavelength,
        'Threshold': header.threshold,
        '# Spots': header.num_spots,
        'Intensity': header.intensity,
        'Ice': header.ice,
        'Max Radius': header.max_distance,
        'SNR': header.snr,
    }

def print_help():
    print("Usage: sporta [options] [files...]")
    print("\nOptions:")
    print("  -h          Display this help message")
    print("  -v          Verbose mode: Save analysis results to 'data.tsv'")
    print("  -r <start>-<end>  Process a range of frames (e.g., 1-100)")
    print("  -i <index>  Set the frame index to process (default: 1)")
    print("  -t <value>  Set the threshold scale factor (default: 6)")
    print("  -e <value>  Set the ring exclusion proximity radius (default: 2)")
    print("  -g <sigma>  Set the Gaussian filtering sigma value")
    print("  -m <size>   Set the minimum pixel size for spots (default: 5)")
    print("  -o          Enable intermediate image output")
    print("\nExamples:")
    print("  sporta -v data.h5")
    print("  sporta -r 1-10 -v data_master.h5")

def main():
    args = sys.argv
    argc = len(args)
    
    # Separate flags and files
    # We need to intercept -r range argument and handle it in Python
    filtered_args = []
    
    range_start = 1
    range_end = 1
    has_range = False
    verbose = False
    
    i = 1
    while i < argc:
        arg = args[i]
        
        if arg == '-h':
            print_help()
            sys.exit(0)
            
        if arg == '-v':
            verbose = True
            i += 1
            continue

        if arg == '-r':
            if i + 1 < argc:
                r_arg = args[i+1]
                try:
                    parts = r_arg.split('-')
                    if len(parts) == 2:
                        range_start = int(parts[0])
                        range_end = int(parts[1])
                        has_range = True
                    else:
                        print(f"Invalid range format: {r_arg}. Expected start-end (e.g. 1-10)")
                        sys.exit(1)
                except ValueError:
                    print(f"Invalid range values: {r_arg}")
                    sys.exit(1)
                i += 2
                continue
            else:
                print("Missing range argument after -r")
                sys.exit(1)
        
        filtered_args.append(arg)
        # Skip values for other flags to avoid double adding
        if arg.startswith('-'):
            # These flags take an argument
            if arg in ['-f', '-i', '-d', '-t', '-s', '-e', '-g', '-m']:
                if i + 1 < argc:
                    filtered_args.append(args[i+1])
                    i += 2
                else:
                    i += 1
            else:
                i += 1
        else:
            i += 1

    # Reconstruct argv for C without -r or -v
    # We will override -i manually if range is set
    c_argv_base = ["sporta"] + filtered_args
    
    # Identify files to process (same logic as before but now using filtered_args)
    files_to_process = []
    i = 0
    while i < len(filtered_args):
        arg = filtered_args[i]
        if arg.startswith('-'):
            if arg in ['-f', '-i', '-d', '-t', '-s', '-e', '-g', '-m']:
                i += 2
            else:
                i += 1
        else:
            files_to_process.append(arg)
            i += 1
            
    # Logic for when files exist OR we want to support range on files
    has_processed = False
    
    if files_to_process:
        collected_results = []
        
        # Determine index list
        indices = range(range_start, range_end + 1) if has_range else [None]
        
        for f in files_to_process:
            if os.path.isdir(f):
                 # Directory logic not fully implemented in this python wrapper for range?
                 # Ignoring for now, focusing on master files
                 pass
            else:
                abs_path = os.path.abspath(f)
                
                for idx in indices:
                    # Construct args for this iteration
                    current_argv = list(c_argv_base)
                    
                    if idx is not None:
                        # Find if -i is already there and replace it, or append it
                        if '-i' in current_argv:
                            idx_pos = current_argv.index('-i')
                            current_argv[idx_pos+1] = str(idx)
                        else:
                            current_argv.extend(['-i', str(idx)])
                            
                    # Update C config
                    c_argc = len(current_argv)
                    c_args = (ctypes.c_char_p * (c_argc + 1))()
                    c_args[:-1] = [arg.encode('utf-8') for arg in current_argv]
                    c_args[-1] = None
                    
                    lib.argparse(c_argc, c_args)
                    
                    res = lib.screener(abs_path.encode())
                    
                    if res == 0: # Success
                        header = lib.getheaderdata()
                        collected_results.append({
                            'File': abs_path,
                            'Frame': idx if idx is not None else 1, # approximate if not set
                            'Resolution': header.resolution,
                            '# Spots': header.num_spots,
                            'Ice': header.ice,
                            'Distance': header.distance,
                            'Wavelength': header.wavelength,
                            'Threshold': header.threshold,
                            'Intensity': header.intensity
                        })
                        has_processed = True

        if collected_results and verbose:
            df = pd.DataFrame(collected_results)
            # Reorder columns slightly for better view
            cols = ['File', 'Frame', 'Resolution', '# Spots', 'Ice', 'Distance', 'Wavelength', 'Threshold', 'Intensity']
            # Only existing columns
            cols = [c for c in cols if c in df.columns]
            print("\nBatch Results:")
            print(df[cols])
            df.to_csv('data.tsv', sep='\t', index=False)
            print("SAVED FILE data.tsv")

    if not has_processed and not files_to_process:
        # Fallback to default directory scanning if no files provided
        if has_range:
             print("Warning: Range -r specified but no files provided. Range ignored for directory scan.")
         
        c_args = (ctypes.c_char_p * (len(c_argv_base) + 1))()
        c_args[:-1] = [arg.encode('utf-8') for arg in c_argv_base]
        c_args[-1] = None
        lib.argparse(len(c_argv_base), c_args)
         
        current_directory = os.getcwd()
        # Use getdirs logic via C or Python? 
        # The previous code had specific logic.
        # But for maintenance, let's rely on C's getrecentfile behavior which iterates directory
        directories = getdirs(current_directory, "data")
        
        data_list = []
        for d in directories:
            header = Hdf5Header()
            lib.getrecentfile(d.encode()) 
            header = lib.getheaderdata()
            
            header_dict = {
                'File': d,
                'Resolution': header.resolution,
                'Distance': header.distance,
                'Wavelength': header.wavelength,
                'Threshold': header.threshold,
                '# Spots': header.num_spots,
                'Intensity': header.intensity,
                'Ice': header.ice,
            }
            data_list.append(header_dict)
            
        if len(data_list) > 0 and verbose:
            df = pd.DataFrame(data_list)
            df.to_csv('data.tsv', sep='\t', index=False)
            print("SAVED FILE data.tsv")

if __name__ == "__main__":
    main()


