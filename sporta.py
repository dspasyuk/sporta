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
        ("postintensity", ctypes.c_double)
    ]

lib.getdirs.restype = ctypes.POINTER(ctypes.c_char_p)
lib.getdirs.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]

lib.getrecentfile.restype = ctypes.c_int
lib.getrecentfile.argtypes = [ctypes.c_char_p]

lib.free_memory.argtypes = [ctypes.POINTER(ctypes.c_char_p), ctypes.c_int]

lib.getheaderdata.restype = Hdf5Header 

def getdirs(currentdir, folder):
    count = ctypes.c_int()
    result = lib.getdirs(currentdir.encode(), folder.encode(), ctypes.byref(count))
    dirs = [result[i].decode() for i in range(count.value)]
    lib.free_memory(result, count)
    
    return dirs

lib.screener.restype = ctypes.c_int
lib.screener.argtypes = [ctypes.c_char_p]

def process_file(filename, index=1, threshold=6, proximity=2, gaussian=None, output=False):
    """
    Process a single HDF5 file or master file.
    
    Args:
        filename (str): Path to the .h5 file.
        index (int): Frame index to process (1-based, default: 1).
        threshold (float): Threshold scale factor (default: 6).
        proximity (float): Ring exclusion proximity radius (default: 2).
        gaussian (float): Gaussian sigma (optional).
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

def main():
    args = sys.argv
    argc = len(args)
    # ... rest of main ...
    
    # Separate flags and files
    # Pass all args to argparse to set flags
    c_args = (ctypes.c_char_p * (argc + 1))()
    c_args[:-1] = [arg.encode('utf-8') for arg in args]
    c_args[-1] = None
    
    lib.argparse(argc, c_args)
    
    # Identify explicit file arguments (ignoring flags and their values)
    # This is a simple heuristic: arguments that don't start with '-' and aren't values of flags
    # But for simplicity, let's just use what argparse didn't consume? 
    # lib.argparse doesn't return consumed args. 
    # Let's check for files in args that exist on disk or match patterns
    
    files_to_process = []
    
    # Skip argv[0] (program name)
    i = 1
    while i < argc:
        arg = args[i]
        if arg.startswith('-'):
            # It's a flag. 
            # Simple assumption: flags like -f, -i, -d, -t, -s, -e, -g take an argument.
            # -h, -o, -r, -p, -n do not.
            if arg in ['-f', '-i', '-d', '-t', '-s', '-e', '-g']:
                i += 2 # Skip flag and value
            else:
                i += 1 # Skip flag
        else:
            # It's potential file or glob
            # Glob expansion is done by shell usually, so sys.argv has the list.
            files_to_process.append(arg)
            i += 1
            
    data_list = []
    
    if files_to_process:
        # Process explicit files
        for f in files_to_process:
            # Check if directory -> use getdirs logic? Or just skip?
            # User intent: 'sporta file.h5' -> process file.
            if os.path.isdir(f):
                # If directory provided as arg, scan it?
                d_files = getdirs(f, "data") # This looks for "data" folder inside or similar?
                # Existing getdirs logic is specific. Let's assume user passes actual files or directories to scan.
                # But getdirs(currentdir, "data") implies folder name to look for.
                pass 
            else:
                # Process file
                # Need to use absolute path for safety if C code expects it?
                abs_path = os.path.abspath(f)
                header = Hdf5Header()
                # screener processes the file and sets header
                # We need to ensure we don't need getrecentfile logic (regex check)
                res = lib.screener(abs_path.encode())
                
                header = lib.getheaderdata()
                header_dict = {
                    'File': abs_path,
                    'Resolution': header.resolution,
                    'Distance': header.distance,
                    'Wavelength': header.wavelength,
                    'Threshold': header.threshold,
                    '# Spots': header.num_spots,
                    'Intensity': header.intensity,
                    'Ice': header.ice,
                }
                data_list.append(header_dict)
                
    else:
        # Fallback to default directory scanning
        currentdir = os.getcwd()
        directories = getdirs(currentdir, "data")
        
        for d in directories:
            header = Hdf5Header()
            lib.getrecentfile(d.encode()) # getrecentfile calls screener internally for 'correct' files
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
            
    df = pd.DataFrame(data_list)
    df = df[df['Resolution'] != 0]
    
    print(df)
    
    namefile = 'data.tsv'
    df.to_csv(namefile, sep='\t', index=False)
    
    print("SAVED FILE", namefile)

if __name__ == "__main__":
    main()


