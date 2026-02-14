from setuptools import setup, Extension
import os
import subprocess
import sys
import glob

# Attempt to find HDF5 via pkg-config
hdf5_include = []
hdf5_lib = []
hdf5_link_args = []

# Function to run pkg-config
def run_pkg_config():
    try:
        pkg_config_cmd = ['pkg-config', '--cflags', '--libs', 'hdf5']
        output = subprocess.check_output(pkg_config_cmd).decode('utf-8').strip()
        for part in output.split():
            if part.startswith('-I'):
                hdf5_include.append(part[2:])
            elif part.startswith('-L'):
                hdf5_lib.append(part[2:])
            elif part.startswith('-l'):
                 pass # setuptools handles libraries
            else:
                 hdf5_link_args.append(part)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

if not run_pkg_config():
    # Fallback to common locations
    if os.path.exists('/usr/include/hdf5/serial'):
        hdf5_include.append('/usr/include/hdf5/serial')
    if os.path.exists('/usr/lib/x86_64-linux-gnu/hdf5/serial'):
         hdf5_lib.append('/usr/lib/x86_64-linux-gnu/hdf5/serial')
    # For HDF5 installed in custom locations (like gsas2main), check env var
    # or fallback to common prefixes
    pass

# Allow explicit override via Environment variable
if os.getenv('HDF5_HOME'):
    h5 = os.getenv('HDF5_HOME')
    hdf5_include.insert(0, os.path.join(h5, 'include'))
    hdf5_lib.insert(0, os.path.join(h5, 'lib'))
    # Also link against hdf5_hl which might be needed
    
module = Extension(
    'sporta_lib',
    sources=[
        'main.c',
        'src/bshuf_h5filter.c',
        'src/bitshuffle.c',
        'src/bitshuffle_core.c',
        'src/iochain.c',
        'src/lz4.c'
    ],
    include_dirs=['src'] + hdf5_include,
    library_dirs=hdf5_lib,
    libraries=['hdf5', 'hdf5_hl', 'm'],
    extra_compile_args=['-Wall', '-fPIC'],
    extra_link_args=hdf5_link_args,
    # We want a shared library that can be loaded by ctypes, 
    # but setuptools builds a module. However, CDLL can load it.
)

setup(
    name='sporta',
    version='1.0',
    description='Spot Detection & Screening of X-Ray Diffraction Images',
    ext_modules=[module],
    py_modules=['sporta'],
    install_requires=['pandas', 'numpy'],
    entry_points={
        'console_scripts': [
            'sporta=sporta:main',
        ],
    },
)
