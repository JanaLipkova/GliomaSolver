/*
 *  MRAGWavelets_StaticData.cpp
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/28/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#include "MRAGCommon.h"
#include "MRAGWavelets_Interp2ndOrder.h"
namespace MRAG
{
const Real Wavelets_Interp2ndOrder::Ha[1] = {1};
const Real Wavelets_Interp2ndOrder::Ga[3] = {-0.5,1.,-0.5};
const Real Wavelets_Interp2ndOrder::Hs[3] = {0.5,1.,0.5};
const Real Wavelets_Interp2ndOrder::Gs[1] = {1};

const int Wavelets_Interp2ndOrder::HaSupport[2] = {0,1};
const int Wavelets_Interp2ndOrder::GaSupport[2] = {-2,1};
const int Wavelets_Interp2ndOrder::HsSupport[2] = {-1,2};
const int Wavelets_Interp2ndOrder::GsSupport[2] = {-1,0};

const Real Wavelets_Interp2ndOrder::PhiSynthesisSupport[2] = {-1.0,1.0};
const Real Wavelets_Interp2ndOrder::PhiAnalysisSupport[2] = {0,1.0};
const Real Wavelets_Interp2ndOrder::CenteringOffset = 0.0;
}

#include "MRAGWavelets_Interp4thOrder.h"
namespace MRAG
{
const Real Wavelets_Interp4thOrder::Ha[1] = {1};
const Real Wavelets_Interp4thOrder::Ga[7] = {1/16.,0,-9/16.,1,-9/16.,0,1/16.};
const Real Wavelets_Interp4thOrder::Hs[7] = {-1/16.,0,9/16.,1,9/16.,0,-1/16.};
const Real Wavelets_Interp4thOrder::Gs[1] = {1};

const int Wavelets_Interp4thOrder::HaSupport[2] = {0,1};
const int Wavelets_Interp4thOrder::GaSupport[2] = {-4,3};
const int Wavelets_Interp4thOrder::HsSupport[2] = {-3,4};
const int Wavelets_Interp4thOrder::GsSupport[2] = {-1,0};

const Real Wavelets_Interp4thOrder::PhiSynthesisSupport[2] = {-2.0,2.0};
const Real Wavelets_Interp4thOrder::PhiAnalysisSupport[2] = {0,1};
const Real Wavelets_Interp4thOrder::CenteringOffset = 0.0;
}

#include "MRAGWavelets_Haar.h"
namespace MRAG
{
const Real Wavelets_Haar::Ha[2] = { +1/2., +1/2. };
const Real Wavelets_Haar::Ga[2] = { +1/2., -1/2. };
const Real Wavelets_Haar::Hs[2] = { +1, +1 };
const Real Wavelets_Haar::Gs[2] = { +1, -1 };

const int Wavelets_Haar::HaSupport[2] = {-1,1};
const int Wavelets_Haar::GaSupport[2] = {-1,1};
const int Wavelets_Haar::HsSupport[2] = {-1,1};
const int Wavelets_Haar::GsSupport[2] = {-1,1};

//const Real Wavelets_Haar::PhiSynthesisSupport[2] = {0.0,2.0};
//const Real Wavelets_Haar::PhiAnalysisSupport[2] = {0,2.};
const Real Wavelets_Haar::CenteringOffset = 0.5;
}

#include "MRAGWavelets_AverageInterp3rdOrder.h"
namespace MRAG
{
	const Real Wavelets_AverageInterp3rdOrder::Ha[2] = { 1/2., 1/2. };
	const Real Wavelets_AverageInterp3rdOrder::Hs[6] = { -1/8.,1/8.,1, 1, 1/8., -1/8.};
	const Real Wavelets_AverageInterp3rdOrder::Ga[6] = { -1/16., -1/16., 1/2., -1/2., 1/16., 1/16. };
	const Real Wavelets_AverageInterp3rdOrder::Gs[2] = { 1, -1};
	
	const int Wavelets_AverageInterp3rdOrder::HaSupport[2] = {-1,1};
	const int Wavelets_AverageInterp3rdOrder::GaSupport[2] = {-3,3};
	const int Wavelets_AverageInterp3rdOrder::HsSupport[2] = {-3,3};
	const int Wavelets_AverageInterp3rdOrder::GsSupport[2] = {-1,1};
	
	const Real Wavelets_AverageInterp3rdOrder::CenteringOffset = 0.5;
}

#include "MRAGWavelets_AverageInterp5thOrder.h"
namespace MRAG
{
	const Real Wavelets_AverageInterp5thOrder::Ha[2] = { 1/2., 1/2. };
	const Real Wavelets_AverageInterp5thOrder::Hs[10] = { 3/128., -3/128., -11./64.,11./64., 1, 1, 11/64., -11./64., -3/128., 3/128.};
	const Real Wavelets_AverageInterp5thOrder::Ga[10] = { 3/256., 3/256., -11/128., -11/128., 1/2., -1/2., 11/128.,11/128.,-3/256., -3/256. };
	const Real Wavelets_AverageInterp5thOrder::Gs[2] = { 1, -1};
	
	const int Wavelets_AverageInterp5thOrder::HaSupport[2] = {-1,1};
	const int Wavelets_AverageInterp5thOrder::GaSupport[2] = {-5,5};
	const int Wavelets_AverageInterp5thOrder::HsSupport[2] = {-5,5};
	const int Wavelets_AverageInterp5thOrder::GsSupport[2] = {-1,1};
	
	const Real Wavelets_AverageInterp5thOrder::CenteringOffset = 0.5;
}