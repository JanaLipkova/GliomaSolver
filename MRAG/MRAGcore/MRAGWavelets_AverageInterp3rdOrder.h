/*
 *  MRAGWavelets_AverageInterp3rdOrder.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 12/9/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

namespace MRAG
{
	struct Wavelets_AverageInterp3rdOrder {
	private:	
		//static 
		static const Real Ha[2];
		static const Real Ga[6];
		static const Real Hs[6];
		static const Real Gs[2];
		
	public:
		static const bool bIsCellCentered = true;
		
		static const int HaSupport[2];
		static const int GaSupport[2];
		static const int HsSupport[2];
		static const int GsSupport[2];
		//static const Real PhiSynthesisSupport[2];
		//static const Real PhiAnalysisSupport[2];
		static const Real CenteringOffset ;
		
		static const Real getCenteringOffset()
		{
			return 0.5;
		}
		
		static const Real getHa(int i) 
		{
			assert(i>=HaSupport[0] && i<HaSupport[1]);
			return Ha[i-HaSupport[0]];
		}
		
		static const Real getGa(int i) 
		{
			assert(i>=GaSupport[0] && i<GaSupport[1]);
			return Ga[i-GaSupport[0]];
		}
		static const Real getHs(int i) 
		{
			assert(i>=HsSupport[0] && i<HsSupport[1]);
			return Hs[i-HsSupport[0]];
		}
		static const Real getGs(int i) 
		{
			assert(i>=GsSupport[0] && i<GsSupport[1]);
			return Gs[i-GsSupport[0]];
		}
	};
}	
