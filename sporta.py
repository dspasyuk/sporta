import ctypes
import sys
import os
import pandas as pd

import glob

# Try to find the setuptools-built library first (sporta_lib*.so)
# This handles cases where it's installed as a package or built in-place
lib_path = None
lib_files = glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), "sporta_lib*.so"))
if lib_files:
    lib_path = lib_files[0]
else:
    # Fallback to manual build name
    manual_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sporta.so")
    if os.path.exists(manual_path):
        lib_path = manual_path

if lib_path is None:
    # If installed via pip, the .so might be in the package directory (which is this directory)
    # but named differently or not found if we are running from source without build
    # Try finding it in site-packages if we are imported?
    # For now, just raise error or assume current directory
     raise FileNotFoundError("Could not find sporta_lib shared library. Did you run 'pip install .' or 'python setup.py build_ext --inplace'?")

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

def main():
    args = sys.argv
    argc = len(args)
    
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


