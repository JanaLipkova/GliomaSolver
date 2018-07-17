
/*
 *  GliomaTypes.h
 *  GliomaXcode
 *
 *  Created by Lipkova on 9/19/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */


#pragma once
#include "Glioma.h"
#include "Matrix.h"

#include "Operators/Glioma_ReactionDiffusionOperator.h"

struct Cell
{
    
    /* tumor */
    Real phi;
    Real dphidt;
    
    /* tissue percentage per voxel*/
	Real p_g, p_w, p_csf;
    
    /* tissue concentration */
    Real wm, gm, csf;  //
    Real dwmdt, dgmdt, dcsfdt;
    
    // velocity field + helper fields for WENO
    Real ux,uy,uz;
    Real omega, domegadt;
    
    // pressure + auxiliary functions for pressure source, phase filed funtion, charact. function
    Real p;
    Real f;
    Real pff, dpffdt;   // pahse field function of whole anatomy, of tissue
    Real chi;
    
    // other helper fields
	Real eps;
	Real tmp;

    
	Cell()
	{
		phi		 = 0.0;
		dphidt	 = 0.0;
		p_g		 = 0.0;
		p_w		 = 0.0;
		p_csf    = 0.0;
        wm = gm = csf = 0.0;
        dwmdt = dgmdt = dcsfdt = 0.0;
		ux = uy = uz  = 0.0;
		p        = 0.0;
		omega	 = 0.0;
		domegadt = 0.0;
		eps      = 0.0;
		tmp      = 0.0;
		f        = 0.0;
		pff		 = 0.0;
		dpffdt   = 0.0;
        chi      = 0.0;

    }
	
    Cell(Real phi_, Real dphidt_, Real p_g_, Real p_w_, Real p_csf_, Real wm_, Real gm_, Real csf_, Real dwmdt_, Real dgmdt_, Real dcsfdt_, Real ux_, Real uy_, Real uz_, Real p_, Real omega_, Real domegadt_, Real eps_, Real tmp_ , Real f_, Real pff_, Real dpffdt_, Real chi_)
	{
		phi		 = phi_	    ;
		dphidt	 = dphidt_  ;
		p_g		 = p_g_	    ;
		p_w		 = p_w_	    ;
		p_csf	 = p_csf_   ;
        wm       = wm_      ;
        gm       = gm_      ;
        csf      = csf_     ;
        dwmdt    = dwmdt_   ;
        dgmdt    = dgmdt_   ;
        dcsfdt   = dcsfdt_  ;
		ux		 = ux_	    ;
		uy		 = uy_	    ;
		uz		 = uz_	    ;
		p		 = p_       ;
		omega    = omega_   ;
		domegadt = domegadt_;
		eps      = eps_     ;
		tmp      = tmp_     ;
		f        = f_       ;
		pff		 = pff_	    ;
		dpffdt   = dpffdt_  ;
        chi      = chi_     ;
    }
	
	void operator += (Cell t)
	{
		phi		 += t.phi	  ;
		dphidt	 += t.dphidt  ;
		p_g		 += t.p_g	  ;
		p_w		 += t.p_w	  ;
		p_csf    += t.p_csf   ;
        wm       += t.wm      ;
        gm       += t.gm      ;
        csf      += t.csf     ;
        dwmdt    += t.dwmdt   ;
        dgmdt    += t.dgmdt   ;
        dcsfdt   += t.dcsfdt  ;
		ux		 += t.ux	  ;
		uy		 += t.uy	  ;
		uz		 += t.uz	  ;
		p		 += t.p 	  ;
		omega    += t.omega	  ;
		domegadt += t.domegadt;
		eps      += t.eps     ;
		tmp      += t.tmp	  ;
		f        += t.f		  ;
		pff		 += t.pff     ;
		dpffdt   += t.dpffdt  ;
        chi      += t.chi     ;
	}
	
	
	
	operator Real() 
	{
		return (Real)phi;
	}
	
	void integrate(float dt)
	{
		phi		+= dt*dphidt;
		dphidt	= 0.0;
	}
	
	template<int i>
	Real evaluate_concentration(double dt)
	{
		return  phi+dt*dphidt;
	}
	
	
	Real giveMe(int i, Real h=0)
	{
		switch(i)
		{
            case 0: return phi;
            case 1: return phi + 0.1 * p_g + 0.2 * p_w;
            case 2: return p_w;
            case 3: return p_g;
            case 4: return p_csf;
            case 5: return chi;
            case 6: return ux;
            case 7: return uy;
            case 8: return uz;

			default: abort(); return 0;
		}
	}
				
	
};

inline Cell operator*(const Cell& p, Real v)
{
	Cell c;
	c.phi		= p.phi		 *v;
	c.dphidt	= p.dphidt	 *v;
	c.p_g		= p.p_g		 *v;
	c.p_w		= p.p_w		 *v;
	c.p_csf     = p.p_csf    *v;
    c.wm        = p.wm       *v;
    c.gm        = p.gm       *v;
    c.csf       = p.csf      *v;
    c.dwmdt     = p.dwmdt    *v;
    c.dgmdt     = p.dgmdt    *v;
    c.dcsfdt    = p.dcsfdt   *v;
	c.ux        = p.ux       *v;
	c.uy        = p.uy       *v;
	c.uz        = p.uz       *v;
	c.p         = p.p        *v;
	c.omega     = p.omega    *v;
	c.domegadt	= p.domegadt *v;
	c.eps       = p.eps      *v;
	c.tmp       = p.tmp      *v;
	c.f         = p.f        *v;
	c.pff		= p.pff		 *v;
	c.dpffdt    = p.dpffdt   *v;
    c.chi       = p.chi      *v;
   
	return c;
}


#pragma mark projectors

template <typename T, int i> 
inline Real RD_projector_impl_vtk(const T&t)
{
	//	return (Real)(0.2*t.p_w + 0.1*t.p_g );	
	return (Real)(0.1 * t.p_g + 0.2 * t.p_w + t.p_csf + t.phi );
}

template <typename T, int i> 
inline Real RD_projector_impl_wav(const T&t)
{
	//return i==0 ? (Real)(t.phi) : (Real)(t.p_w);  // for refinment w.r.t 2 channels
#ifdef LiverModel
    return (Real)(t.psi) ;
#else
    return (Real)(t.phi) ;
#endif

}

make_projector(RD_Projector_VTK,      RD_projector_impl_vtk)
make_projector(RD_Projector_Wavelets, RD_projector_impl_wav)


#ifndef _DIM
#define _DIM 3
#endif

#ifndef _BLOCKSIZE_
#define _BLOCKSIZE_ 16
#endif

#ifndef _BLOCKSIZE_Z_
#define _BLOCKSIZE_Z_ _BLOCKSIZE_
#endif

#ifndef _BPD_
#define _BPD_ 4
#endif

#ifndef _MAXLEVEL_
#define _MAXLEVEL_ 3
#endif

static const int blockSize = _BLOCKSIZE_;
static const int blockSizeZ = _BLOCKSIZE_Z_;
static const int blocksPerDimension = _BPD_;


// Structural parameters
static const bool	bIsCellCentered = true;
static const bool	bVerbose		= true;

// Multiresolution parameters
static const int maxLevel = _MAXLEVEL_;    // 8bpd has maxLevel 3 since 2^3
static const int resJump  = 1;    // modulo(maxLevel,resJum) = 0, !!! and reJump < maxLevel
const double refinement_tolerance	= 1e-3;
const double compression_tolerance	= 1e-5;

typedef		Block< Cell, blockSize, blockSize, blockSizeZ>					B;
typedef		_WAVELET_TYPE															W;

typedef Matrix::D3D<double> MatrixD3D;
typedef Matrix::D2D<double> MatrixD2D;

#ifdef _MRAG_TBB
static const int nThreads = _MRAG_TBB_NTHREADS_HINT ;
typedef		Multithreading::BlockProcessing_Pipeline_TBB<B, BlockLab, nThreads+1>	BlockProcessing;
#else
static const int nThreads = 1 ;
typedef		Multithreading::BlockProcessing_SingleCPU<B>							BlockProcessing;
#endif
