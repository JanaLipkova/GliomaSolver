LIB_BASE="_USER_LIB_BASE_"

export DYLD_LIBRARY_PATH=${LIB_BASE}/tbb40_20120613oss/build/macos_intel64_gcc_cc4.8.5_os10.10.5_release/:$DYLD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=${LIB_BASE}/hypre-2.10.0b/src/hypre/lib/:$DYLD_LIBRARY_PATH

# For macports installation:
export DYLD_LIBRARY_PATH=/opt/local/lib/vtk-5.10/:$DYLD_LIBRARY_PATH
# Default location if installed from source:
# export DYLD_LIBRARY_PATH=/usr/local/lib/vtk-5.6/:$DYLD_LIBRARY_PATH

