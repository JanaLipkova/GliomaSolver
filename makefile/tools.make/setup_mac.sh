LIB_BASE="_USER_LIB_BASE_"

export DYLD_LIBRARY_PATH=${LIB_BASE}/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:$DYLD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=${LIB_BASE}/myVTK/lib/vtk-5.4/:$DYLD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=${LIB_BASE}/hypre-2.10.0b/src/hypre/lib/:$DYLD_LIBRARY_PATH

