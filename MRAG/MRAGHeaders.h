/*
 *  MRAGHeaders.h
 *  GliomaCode
 *
 *  Created by Diego Rossinelli on 7/15/10, 
 *	modified by Jana Lipkova 
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGEnvironment.h"

#define _WAVELET_TYPE Wavelets_Interp2ndOrder
#undef _MRAG_GLUT_VIZ


#if _MRAG_OS == _MRAG_OS_APPLE
#ifdef _MRAG_GLUT_VIZ
#include "GLUT/glut.h"
#endif
#elif _MRAG_OS == _MRAG_OS_WINDOWS
#define _USE_MATH_DEFINES
#ifdef _MRAG_GLUT_VIZ
#include "GL/glew.h"
#include "GL/glut.h"
#endif
#endif

#undef min
#undef max


#include "MRAGcore/MRAGWavelets_AverageInterp5thOrder.h"
#include "MRAGcore/MRAGWavelets_Interp4thOrder.h"
#include "MRAGcore/MRAGWavelets_AverageInterp3rdOrder.h"
#include "MRAGcore/MRAGWavelets_Interp2ndOrder.h"
#include "MRAGcore/MRAGWavelets_Haar.h"
#include "MRAGcore/MRAGrid.h"
#include "MRAGcore/MRAGRefiner.h"
#include "MRAGcore/MRAGCompressor.h"
#include "MRAGcore/MRAGBlockLab.h"
#include "MRAGcore/MRAGBlockFWT.h"
#include "MRAGcore/MRAGProfiler.h"

#ifdef _MRAG_GLUT_VIZ
    #include "MRAGvisual/GridViewer.h"
#endif

#include "MRAGscience/MRAGScienceCore.h"
#include "MRAGscience/MRAGAutomaticRefiner.h"
#include "MRAGscience/MRAGSimpleLevelsetBlock.h"
#include "MRAGscience/MRAGSpaceTimeSorter.h"
#include "MRAGscience/MRAGRefiner_SpaceExtension.h"
#include "MRAGscience/MRAGWeno.h"

#include "MRAGmultithreading/MRAGBlockProcessing_SingleCPU.h"
#ifdef _MRAG_TBB
    #include "MRAGmultithreading/MRAGBlockProcessing_TBB.h"
#endif

#include "MRAGio/MRAG_IO_ArgumentParser.h"
#include "MRAGio/MRAG_IO_Native.h"
#include "MRAGio/MRAG_IO_VTK.h"
#include "MRAGio/MRAG_IO_VTKNative.h"
#include "MRAGio/MRAG_IO_VTU3D.h"
//#include "MRAGio/MRAG_IO_VTKMB.h"
//#include "MRAGio/DumpScalarToVP.h"

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <omp.h>
#include <stdlib.h>     /* abs, rand */
#include <math.h>       /* fabs */
