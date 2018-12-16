LIB_BASE="_USER_LIB_BASE_"

export DYLD_LIBRARY_PATH=${LIB_BASE}/tbb40_20120613oss/build/macos_intel64_gcc_cc4.8.5_os10.10.5_release/:$DYLD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=${LIB_BASE}/hypre-2.10.0b/src/hypre/lib/:$DYLD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=/opt/local/lib/vtk-5.10/:$DYLD_LIBRARY_PATH

