#pragma once
#include <math.h>
#include <assert.h>

#include "MRAGCommon.h"


namespace MRAG
{
class Wavelets_Bior22 {
private:
	static Real _Kernel(Real x)
	{
		return max(0.0, 1.0-fabs(x));
	}

	static const Real Ha[5];
	static const Real Ga[3];
	static const Real Hs[3];
	static const Real Gs[5];

public:
	static const int HaSupport[2];
	static const int GaSupport[2];
	static const int HsSupport[2];
	static const int GsSupport[2];
	static const Real PhiSynthesisSupport[2]; 
	
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

	static const Real PhiSynthesis(Real x) 
	{
		return fabs(x)<0.5?1:0;
		//return _Kernel(x);
	}
	
	static const Real PhiAnalysis(Real x) 
	{
		return fabs(x)<0.5?1:0;
		//return _Kernel(x);
	}
	static const Real PhiSynthesis(int level, int offset, Real x) 
	{
		return pow(2.0, level*0.5)*_Kernel(pow(2.0, level)*x - offset);
	}

	static const Real PsiSynthesis(int level, int offset, Real x) 
	{
		Real val = 0;

		for(int i=GsSupport[0]; i<GsSupport[1]; i++)
			val += getGs(i)*PhiSynthesis(level+1, i + 2*offset, x);

		return val;
	}
};	


const Real Wavelets_Bior22::Ha[5] = {-0.176776695296637,   0.353553390593274,   1.060660171779821,   0.353553390593274,  -0.176776695296637};
const Real Wavelets_Bior22::Ga[3] = {0.353553390593274,  -0.707106781186548,   0.353553390593274};
const Real Wavelets_Bior22::Hs[3] = {0.353553390593274,   0.707106781186548,   0.353553390593274};
const Real Wavelets_Bior22::Gs[5] = {0.176776695296637,   0.353553390593274,  -1.060660171779821,   0.353553390593274,   0.176776695296637};

const int Wavelets_Bior22::HaSupport[2] = {-2,3};
const int Wavelets_Bior22::GaSupport[2] = {0,3};
const int Wavelets_Bior22::HsSupport[2] = {-1,2};
const int Wavelets_Bior22::GsSupport[2] = {-1,4};

const Real Wavelets_Bior22::PhiSynthesisSupport[2] = {-0.5,0.5}; 
}

	
