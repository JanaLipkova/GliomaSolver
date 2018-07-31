//
//  Glioma_TissueAdvectionOperator.h
//  GliomaXcode
//
//  Created by Lipkova on 21/05/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#ifndef GliomaXcode_Glioma_TissueAdvectionOperator_h
#define GliomaXcode_Glioma_TissueAdvectionOperator_h


#define NoFluxBC
#define CASE1

template<int nDim = 3>
struct TissueAdvectionWeno5
{
    int stencil_start[3];
    int stencil_end[3];
    
    TissueAdvectionWeno5()
    {
        stencil_start[0]=stencil_start[1] = -3;
        stencil_end[0]  =stencil_end[1]   = +4;
        stencil_start[2] = nDim==3 ? -3: 0;
        stencil_end[2]   = nDim==3 ? +4:+1;
    }
    
    TissueAdvectionWeno5(const TissueAdvectionWeno5&)
    {
        stencil_start[0]=stencil_start[1] = -3;
        stencil_end[0]  =stencil_end[1]   = +4;
        stencil_start[2] = nDim==3 ? -3: 0;
        stencil_end[2]   = nDim==3 ? +4:+1;
        
    }
    
    template<typename LabType, typename BlockType>
    inline void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        Weno5<double> weno5;
        
        const double ifactor = 1./(info.h[0]);
        
        // fluxes
        double wm_flux_x[2];
        double wm_flux_y[2];
        double wm_flux_z[2];
        
        double gm_flux_x[2];
        double gm_flux_y[2];
        double gm_flux_z[2];
        
        double csf_flux_x[2];
        double csf_flux_y[2];
        double csf_flux_z[2];
        
        double tumor_flux_x[2];
        double tumor_flux_y[2];
        double tumor_flux_z[2];
        
        
        Real n[6]; // neighbours [-3,-2,-1,1,2,3]
        
        if(nDim==2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    if( lab(ix,iy).chi > 0)
                    {
                        /* In X */
                        n[0] = lab(ix-3,iy).chi;
                        n[1] = lab(ix-2,iy).chi;
                        n[2] = lab(ix-1,iy).chi;
                        n[3] = lab(ix+1,iy).chi;
                        n[4] = lab(ix+2,iy).chi;
                        n[5] = lab(ix+3,iy).chi;
                        
                        // _applyNoFluxBC(n);
                        
                        // white matter
                        weno5(  ifactor * ( n[1] * lab(ix-2,iy).wm - n[0] * lab(ix-3,iy).wm ),
                              ifactor   * ( n[2] * lab(ix-1,iy).wm - n[1] * lab(ix-2,iy).wm ),
                              ifactor   * (        lab(ix  ,iy).wm - n[2] * lab(ix-1,iy).wm ),
                              ifactor   * ( n[3] * lab(ix+1,iy).wm -        lab(ix  ,iy).wm ),
                              ifactor   * ( n[4] * lab(ix+2,iy).wm - n[3] * lab(ix+1,iy).wm ),
                              ifactor   * ( n[5] * lab(ix+3,iy).wm - n[4] * lab(ix+2,iy).wm ),
                              wm_flux_x );
                        
                        // gray matter
                        weno5(  ifactor * ( n[1] * lab(ix-2,iy).gm - n[0] * lab(ix-3,iy).gm ),
                              ifactor   * ( n[2] * lab(ix-1,iy).gm - n[1] * lab(ix-2,iy).gm ),
                              ifactor   * (        lab(ix  ,iy).gm - n[2] * lab(ix-1,iy).gm ),
                              ifactor   * ( n[3] * lab(ix+1,iy).gm -        lab(ix  ,iy).gm ),
                              ifactor   * ( n[4] * lab(ix+2,iy).gm - n[3] * lab(ix+1,iy).gm ),
                              ifactor   * ( n[5] * lab(ix+3,iy).gm - n[4] * lab(ix+2,iy).gm ),
                              gm_flux_x );
                        
                        // csf
                        weno5(  ifactor * ( n[1] * lab(ix-2,iy).csf - n[0] * lab(ix-3,iy).csf ),
                              ifactor   * ( n[2] * lab(ix-1,iy).csf - n[1] * lab(ix-2,iy).csf ),
                              ifactor   * (        lab(ix  ,iy).csf - n[2] * lab(ix-1,iy).csf ),
                              ifactor   * ( n[3] * lab(ix+1,iy).csf -        lab(ix  ,iy).csf ),
                              ifactor   * ( n[4] * lab(ix+2,iy).csf - n[3] * lab(ix+1,iy).csf ),
                              ifactor   * ( n[5] * lab(ix+3,iy).csf - n[4] * lab(ix+2,iy).csf ),
                              csf_flux_x );
                        
                        // tumor
                        n[0] = lab(ix-3,iy).p_w + lab(ix-3,iy).p_g ;
                        n[1] = lab(ix-2,iy).p_w + lab(ix-2,iy).p_g ;
                        n[2] = lab(ix-1,iy).p_w + lab(ix-1,iy).p_g ;
                        n[3] = lab(ix+1,iy).p_w + lab(ix+1,iy).p_g ;
                        n[4] = lab(ix+2,iy).p_w + lab(ix+2,iy).p_g ;
                        n[5] = lab(ix+3,iy).p_w + lab(ix+3,iy).p_g ;

                        
#ifdef NoFluxBC
                        _applyNoFluxBC(n);
#endif
                        
                        weno5(  ifactor * ( n[1] * lab(ix-2,iy).phi - n[0] * lab(ix-3,iy).phi ),
                              ifactor   * ( n[2] * lab(ix-1,iy).phi - n[1] * lab(ix-2,iy).phi ),
                              ifactor   * (        lab(ix  ,iy).phi - n[2] * lab(ix-1,iy).phi ),
                              ifactor   * ( n[3] * lab(ix+1,iy).phi -        lab(ix  ,iy).phi ),
                              ifactor   * ( n[4] * lab(ix+2,iy).phi - n[3] * lab(ix+1,iy).phi ),
                              ifactor   * ( n[5] * lab(ix+3,iy).phi - n[4] * lab(ix+2,iy).phi ),
                              tumor_flux_x );
                        
                        
                        /* In Y */
                        n[0] = lab(ix,iy-3).chi;
                        n[1] = lab(ix,iy-2).chi;
                        n[2] = lab(ix,iy-1).chi;
                        n[3] = lab(ix,iy+1).chi;
                        n[4] = lab(ix,iy+2).chi;
                        n[5] = lab(ix,iy+3).chi;
                        
                        //_applyNoFluxBC(n);
                        
                        // white matter
                        weno5( ifactor * ( n[1] * lab(ix,iy-2).wm - n[0] * lab(ix,iy-3).wm ),
                              ifactor  * ( n[2] * lab(ix,iy-1).wm - n[1] * lab(ix,iy-2).wm ),
                              ifactor  * (        lab(ix,iy  ).wm - n[2] * lab(ix,iy-1).wm ),
                              ifactor  * ( n[3] * lab(ix,iy+1).wm -        lab(ix,iy  ).wm ),
                              ifactor  * ( n[4] * lab(ix,iy+2).wm - n[3] * lab(ix,iy+1).wm ),
                              ifactor  * ( n[5] * lab(ix,iy+3).wm - n[4] * lab(ix,iy+2).wm ),
                              wm_flux_y );
                        
                        // gray matter
                        weno5( ifactor * ( n[1] * lab(ix,iy-2).gm - n[0] * lab(ix,iy-3).gm ),
                              ifactor  * ( n[2] * lab(ix,iy-1).gm - n[1] * lab(ix,iy-2).gm ),
                              ifactor  * (        lab(ix,iy  ).gm - n[2] * lab(ix,iy-1).gm ),
                              ifactor  * ( n[3] * lab(ix,iy+1).gm -        lab(ix,iy  ).gm ),
                              ifactor  * ( n[4] * lab(ix,iy+2).gm - n[3] * lab(ix,iy+1).gm ),
                              ifactor  * ( n[5] * lab(ix,iy+3).gm - n[4] * lab(ix,iy+2).gm ),
                              gm_flux_y );
                        
                        // csf
                        weno5( ifactor * ( n[1] * lab(ix,iy-2).csf - n[0] * lab(ix,iy-3).csf ),
                              ifactor  * ( n[2] * lab(ix,iy-1).csf - n[1] * lab(ix,iy-2).csf ),
                              ifactor  * (        lab(ix,iy  ).csf - n[2] * lab(ix,iy-1).csf ),
                              ifactor  * ( n[3] * lab(ix,iy+1).csf -        lab(ix,iy  ).csf ),
                              ifactor  * ( n[4] * lab(ix,iy+2).csf - n[3] * lab(ix,iy+1).csf ),
                              ifactor  * ( n[5] * lab(ix,iy+3).csf - n[4] * lab(ix,iy+2).csf ),
                              csf_flux_y );
                        
                        // tumor
                        n[0] = lab(ix,iy-3).p_w + lab(ix,iy-3).p_g ;
                        n[1] = lab(ix,iy-2).p_w + lab(ix,iy-2).p_g ;
                        n[2] = lab(ix,iy-1).p_w + lab(ix,iy-1).p_g ;
                        n[3] = lab(ix,iy+1).p_w + lab(ix,iy+1).p_g ;
                        n[4] = lab(ix,iy+2).p_w + lab(ix,iy+2).p_g ;
                        n[5] = lab(ix,iy+3).p_w + lab(ix,iy+3).p_g ;
                        
#ifdef NoFluxBC
                        _applyNoFluxBC(n);
#endif
                        
                        
                        
                        weno5( ifactor * ( n[1] * lab(ix,iy-2).phi - n[0] * lab(ix,iy-3).phi ),
                              ifactor  * ( n[2] * lab(ix,iy-1).phi - n[1] * lab(ix,iy-2).phi ),
                              ifactor  * (        lab(ix,iy  ).phi - n[2] * lab(ix,iy-1).phi ),
                              ifactor  * ( n[3] * lab(ix,iy+1).phi -        lab(ix,iy  ).phi ),
                              ifactor  * ( n[4] * lab(ix,iy+2).phi - n[3] * lab(ix,iy+1).phi ),
                              ifactor  * ( n[5] * lab(ix,iy+3).phi - n[4] * lab(ix,iy+2).phi ),
                              tumor_flux_y );
                        
                        // velocities
                        Real velocity_plus[2] = {
                            max(o(ix,iy).ux, (Real) 0.),
                            max(o(ix,iy).uy, (Real) 0.)
                        };
                        
                        Real velocity_minus[2] = {
                            min(lab(ix,iy).ux, (Real) 0.),
                            min(lab(ix,iy).uy, (Real) 0.)
                        };
                        
                        o(ix,iy).dwmdt  = (-1.) * ( velocity_plus[0] * wm_flux_x[0] + velocity_minus[0] * wm_flux_x[1] +
                                                   velocity_plus[1] * wm_flux_y[0] + velocity_minus[1] * wm_flux_y[1]  );
                        
                        o(ix,iy).dgmdt  = (-1.) * ( velocity_plus[0] * gm_flux_x[0] + velocity_minus[0] * gm_flux_x[1] +
                                                   velocity_plus[1] * gm_flux_y[0] + velocity_minus[1] * gm_flux_y[1]  );
                        
                        o(ix,iy).dcsfdt = (-1.) * ( velocity_plus[0] * csf_flux_x[0] + velocity_minus[0] * csf_flux_x[1] +
                                                   velocity_plus[1] * csf_flux_y[0] + velocity_minus[1] * csf_flux_y[1]  );
                        
                        o(ix,iy).dphidt = (-1.) * ( velocity_plus[0] * tumor_flux_x[0] + velocity_minus[0] * tumor_flux_x[1] +
                                                   velocity_plus[1] * tumor_flux_y[0] + velocity_minus[1] * tumor_flux_y[1]  );
                        
                    }
                    else
                    {
                        o(ix,iy).dwmdt  = 0.;
                        o(ix,iy).dgmdt  = 0.;
                        o(ix,iy).dcsfdt = 0.;
                        o(ix,iy).dphidt = 0.;
                        
                    }
                    
                }
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        if( lab(ix,iy,iz).chi > 0)
                        {
                            /* In X */
                            n[0] = lab(ix-3,iy,iz).chi;
                            n[1] = lab(ix-2,iy,iz).chi;
                            n[2] = lab(ix-1,iy,iz).chi;
                            n[3] = lab(ix+1,iy,iz).chi;
                            n[4] = lab(ix+2,iy,iz).chi;
                            n[5] = lab(ix+3,iy,iz).chi;
                            
                            // _applyNoFluxBC(n);
                            
                            // white matter
                            weno5(  ifactor * ( n[1] * lab(ix-2,iy,iz).wm - n[0] * lab(ix-3,iy,iz).wm ),
                                  ifactor   * ( n[2] * lab(ix-1,iy,iz).wm - n[1] * lab(ix-2,iy,iz).wm ),
                                  ifactor   * (        lab(ix  ,iy,iz).wm - n[2] * lab(ix-1,iy,iz).wm ),
                                  ifactor   * ( n[3] * lab(ix+1,iy,iz).wm -        lab(ix  ,iy,iz).wm ),
                                  ifactor   * ( n[4] * lab(ix+2,iy,iz).wm - n[3] * lab(ix+1,iy,iz).wm ),
                                  ifactor   * ( n[5] * lab(ix+3,iy,iz).wm - n[4] * lab(ix+2,iy,iz).wm ),
                                  wm_flux_x );
                            
                            // gray matter
                            weno5(  ifactor * ( n[1] * lab(ix-2,iy,iz).gm - n[0] * lab(ix-3,iy,iz).gm ),
                                  ifactor   * ( n[2] * lab(ix-1,iy,iz).gm - n[1] * lab(ix-2,iy,iz).gm ),
                                  ifactor   * (        lab(ix  ,iy,iz).gm - n[2] * lab(ix-1,iy,iz).gm ),
                                  ifactor   * ( n[3] * lab(ix+1,iy,iz).gm -        lab(ix  ,iy,iz).gm ),
                                  ifactor   * ( n[4] * lab(ix+2,iy,iz).gm - n[3] * lab(ix+1,iy,iz).gm ),
                                  ifactor   * ( n[5] * lab(ix+3,iy,iz).gm - n[4] * lab(ix+2,iy,iz).gm ),
                                  gm_flux_x );
                            
                            // csf
                            weno5(  ifactor * ( n[1] * lab(ix-2,iy,iz).csf - n[0] * lab(ix-3,iy,iz).csf ),
                                  ifactor   * ( n[2] * lab(ix-1,iy,iz).csf - n[1] * lab(ix-2,iy,iz).csf ),
                                  ifactor   * (        lab(ix  ,iy,iz).csf - n[2] * lab(ix-1,iy,iz).csf ),
                                  ifactor   * ( n[3] * lab(ix+1,iy,iz).csf -        lab(ix  ,iy,iz).csf ),
                                  ifactor   * ( n[4] * lab(ix+2,iy,iz).csf - n[3] * lab(ix+1,iy,iz).csf ),
                                  ifactor   * ( n[5] * lab(ix+3,iy,iz).csf - n[4] * lab(ix+2,iy,iz).csf ),
                                  csf_flux_x );
                            
                           
#ifdef TissueMemory
                            /* In X for Tumour*/
                            n[0] = lab(ix-3,iy,iz).p_w + lab(ix-3,iy,iz).p_g ;
                            n[1] = lab(ix-2,iy,iz).p_w + lab(ix-2,iy,iz).p_g ;
                            n[2] = lab(ix-1,iy,iz).p_w + lab(ix-1,iy,iz).p_g ;
                            n[3] = lab(ix+1,iy,iz).p_w + lab(ix+1,iy,iz).p_g ;
                            n[4] = lab(ix+2,iy,iz).p_w + lab(ix+2,iy,iz).p_g ;
                            n[5] = lab(ix+3,iy,iz).p_w + lab(ix+3,iy,iz).p_g ;
#else
      
                            n[0] = lab(ix-3,iy,iz).p_w + lab(ix-3,iy,iz).p_g + lab(ix-3,iy,iz).phi;
                            n[1] = lab(ix-2,iy,iz).p_w + lab(ix-2,iy,iz).p_g + lab(ix-2,iy,iz).phi;
                            n[2] = lab(ix-1,iy,iz).p_w + lab(ix-1,iy,iz).p_g + lab(ix-1,iy,iz).phi;
                            n[3] = lab(ix+1,iy,iz).p_w + lab(ix+1,iy,iz).p_g + lab(ix+1,iy,iz).phi;
                            n[4] = lab(ix+2,iy,iz).p_w + lab(ix+2,iy,iz).p_g + lab(ix+2,iy,iz).phi;
                            n[5] = lab(ix+3,iy,iz).p_w + lab(ix+3,iy,iz).p_g + lab(ix+3,iy,iz).phi;
                            
                            // normalise
                            n[0] = n[0] / (n[0] + lab(ix-3,iy,iz).p_csf );
                            n[1] = n[1] / (n[1] + lab(ix-2,iy,iz).p_csf );
                            n[2] = n[2] / (n[2] + lab(ix-1,iy,iz).p_csf );
                            n[3] = n[3] / (n[3] + lab(ix+1,iy,iz).p_csf );
                            n[4] = n[4] / (n[4] + lab(ix+2,iy,iz).p_csf );
                            n[5] = n[5] / (n[5] + lab(ix+3,iy,iz).p_csf );
#endif
                            
#ifdef NoFluxBC
                            _applyNoFluxBC(n);
#endif
                            
                            // tumor
                            weno5(  ifactor * ( n[1] * lab(ix-2,iy,iz).phi - n[0] * lab(ix-3,iy,iz).phi ),
                                  ifactor   * ( n[2] * lab(ix-1,iy,iz).phi - n[1] * lab(ix-2,iy,iz).phi ),
                                  ifactor   * (        lab(ix  ,iy,iz).phi - n[2] * lab(ix-1,iy,iz).phi ),
                                  ifactor   * ( n[3] * lab(ix+1,iy,iz).phi -        lab(ix  ,iy,iz).phi ),
                                  ifactor   * ( n[4] * lab(ix+2,iy,iz).phi - n[3] * lab(ix+1,iy,iz).phi ),
                                  ifactor   * ( n[5] * lab(ix+3,iy,iz).phi - n[4] * lab(ix+2,iy,iz).phi ),
                                  tumor_flux_x );
                            
                            
                            
                            
                            /* In Y */
                            n[0] = lab(ix,iy-3,iz).chi;
                            n[1] = lab(ix,iy-2,iz).chi;
                            n[2] = lab(ix,iy-1,iz).chi;
                            n[3] = lab(ix,iy+1,iz).chi;
                            n[4] = lab(ix,iy+2,iz).chi;
                            n[5] = lab(ix,iy+3,iz).chi;
                            
                            // _applyNoFluxBC(n);
                            
                            // white matter
                            weno5( ifactor * ( n[1] * lab(ix,iy-2,iz).wm - n[0] * lab(ix,iy-3,iz).wm ),
                                  ifactor  * ( n[2] * lab(ix,iy-1,iz).wm - n[1] * lab(ix,iy-2,iz).wm ),
                                  ifactor  * (        lab(ix,iy  ,iz).wm - n[2] * lab(ix,iy-1,iz).wm ),
                                  ifactor  * ( n[3] * lab(ix,iy+1,iz).wm -        lab(ix,iy  ,iz).wm ),
                                  ifactor  * ( n[4] * lab(ix,iy+2,iz).wm - n[3] * lab(ix,iy+1,iz).wm ),
                                  ifactor  * ( n[5] * lab(ix,iy+3,iz).wm - n[4] * lab(ix,iy+2,iz).wm ),
                                  wm_flux_y );
                            
                            // gray matter
                            weno5( ifactor * ( n[1] * lab(ix,iy-2,iz).gm - n[0] * lab(ix,iy-3,iz).gm ),
                                  ifactor  * ( n[2] * lab(ix,iy-1,iz).gm - n[1] * lab(ix,iy-2,iz).gm ),
                                  ifactor  * (        lab(ix,iy  ,iz).gm - n[2] * lab(ix,iy-1,iz).gm ),
                                  ifactor  * ( n[3] * lab(ix,iy+1,iz).gm -        lab(ix,iy  ,iz).gm ),
                                  ifactor  * ( n[4] * lab(ix,iy+2,iz).gm - n[3] * lab(ix,iy+1,iz).gm ),
                                  ifactor  * ( n[5] * lab(ix,iy+3,iz).gm - n[4] * lab(ix,iy+2,iz).gm ),
                                  gm_flux_y );
                            
                            // csf
                            weno5( ifactor * ( n[1] * lab(ix,iy-2,iz).csf - n[0] * lab(ix,iy-3,iz).csf ),
                                  ifactor  * ( n[2] * lab(ix,iy-1,iz).csf - n[1] * lab(ix,iy-2,iz).csf ),
                                  ifactor  * (        lab(ix,iy  ,iz).csf - n[2] * lab(ix,iy-1,iz).csf ),
                                  ifactor  * ( n[3] * lab(ix,iy+1,iz).csf -        lab(ix,iy  ,iz).csf ),
                                  ifactor  * ( n[4] * lab(ix,iy+2,iz).csf - n[3] * lab(ix,iy+1,iz).csf ),
                                  ifactor  * ( n[5] * lab(ix,iy+3,iz).csf - n[4] * lab(ix,iy+2,iz).csf ),
                                  csf_flux_y );
                           
#ifdef TissueMemory
                            /* In Y for Tumour*/
                            n[0] = lab(ix,iy-3,iz).p_w + lab(ix,iy-3,iz).p_g ;
                            n[1] = lab(ix,iy-2,iz).p_w + lab(ix,iy-2,iz).p_g ;
                            n[2] = lab(ix,iy-1,iz).p_w + lab(ix,iy-1,iz).p_g ;
                            n[3] = lab(ix,iy+1,iz).p_w + lab(ix,iy+1,iz).p_g ;
                            n[4] = lab(ix,iy+2,iz).p_w + lab(ix,iy+2,iz).p_g ;
                            n[5] = lab(ix,iy+3,iz).p_w + lab(ix,iy+3,iz).p_g ;
#else
      
                            n[0] = lab(ix,iy-3,iz).p_w + lab(ix,iy-3,iz).p_g + lab(ix,iy-3,iz).phi ;
                            n[1] = lab(ix,iy-2,iz).p_w + lab(ix,iy-2,iz).p_g + lab(ix,iy-2,iz).phi ;
                            n[2] = lab(ix,iy-1,iz).p_w + lab(ix,iy-1,iz).p_g + lab(ix,iy-1,iz).phi ;
                            n[3] = lab(ix,iy+1,iz).p_w + lab(ix,iy+1,iz).p_g + lab(ix,iy+1,iz).phi ;
                            n[4] = lab(ix,iy+2,iz).p_w + lab(ix,iy+2,iz).p_g + lab(ix,iy+2,iz).phi ;
                            n[5] = lab(ix,iy+3,iz).p_w + lab(ix,iy+3,iz).p_g + lab(ix,iy+3,iz).phi ;
                            
                            // normalised
                            n[0] = n[0] / ( n[0] + lab(ix,iy-3,iz).p_csf );
                            n[1] = n[1] / ( n[1] + lab(ix,iy-2,iz).p_csf );
                            n[2] = n[2] / ( n[2] + lab(ix,iy-1,iz).p_csf );
                            n[3] = n[3] / ( n[3] + lab(ix,iy+1,iz).p_csf );
                            n[4] = n[4] / ( n[4] + lab(ix,iy+2,iz).p_csf );
                            n[5] = n[5] / ( n[5] + lab(ix,iy+3,iz).p_csf );
#endif
                            
                            

#ifdef NoFluxBC
                            _applyNoFluxBC(n);
#endif
                            
                            
                            // tumor
                            weno5( ifactor * ( n[1] * lab(ix,iy-2,iz).phi - n[0] * lab(ix,iy-3,iz).phi ),
                                  ifactor  * ( n[2] * lab(ix,iy-1,iz).phi - n[1] * lab(ix,iy-2,iz).phi ),
                                  ifactor  * (        lab(ix,iy  ,iz).phi - n[2] * lab(ix,iy-1,iz).phi ),
                                  ifactor  * ( n[3] * lab(ix,iy+1,iz).phi -        lab(ix,iy  ,iz).phi ),
                                  ifactor  * ( n[4] * lab(ix,iy+2,iz).phi - n[3] * lab(ix,iy+1,iz).phi ),
                                  ifactor  * ( n[5] * lab(ix,iy+3,iz).phi - n[4] * lab(ix,iy+2,iz).phi ),
                                  tumor_flux_y );
                            
                            
                            
                            /* In Z */
                            n[0] = lab(ix,iy,iz-3).chi;
                            n[1] = lab(ix,iy,iz-2).chi;
                            n[2] = lab(ix,iy,iz-1).chi;
                            n[3] = lab(ix,iy,iz+1).chi;
                            n[4] = lab(ix,iy,iz+2).chi;
                            n[5] = lab(ix,iy,iz+3).chi;
                            
                            // _applyNoFluxBC(n);
                            
                            // white matter
                            weno5( ifactor * ( n[1] * lab(ix,iy,iz-2).wm - n[0] * lab(ix,iy,iz-3).wm ),
                                  ifactor  * ( n[2] * lab(ix,iy,iz-1).wm - n[1] * lab(ix,iy,iz-2).wm ),
                                  ifactor  * (        lab(ix,iy,iz  ).wm - n[2] * lab(ix,iy,iz-1).wm ),
                                  ifactor  * ( n[3] * lab(ix,iy,iz+1).wm -        lab(ix,iy,iz  ).wm ),
                                  ifactor  * ( n[4] * lab(ix,iy,iz+2).wm - n[3] * lab(ix,iy,iz+1).wm ),
                                  ifactor  * ( n[5] * lab(ix,iy,iz+3).wm - n[4] * lab(ix,iy,iz+2).wm ),
                                  wm_flux_z );
                            
                            // gray matter
                            weno5( ifactor * ( n[1] * lab(ix,iy,iz-2).gm - n[0] * lab(ix,iy,iz-3).gm ),
                                  ifactor  * ( n[2] * lab(ix,iy,iz-1).gm - n[1] * lab(ix,iy,iz-2).gm ),
                                  ifactor  * (        lab(ix,iy  ,iz).gm - n[2] * lab(ix,iy,iz-1).gm ),
                                  ifactor  * ( n[3] * lab(ix,iy,iz+1).gm -        lab(ix,iy,iz  ).gm ),
                                  ifactor  * ( n[4] * lab(ix,iy,iz+2).gm - n[3] * lab(ix,iy,iz+1).gm ),
                                  ifactor  * ( n[5] * lab(ix,iy,iz+3).gm - n[4] * lab(ix,iy,iz+2).gm ),
                                  gm_flux_z );
                            
                            // csf
                            weno5( ifactor * ( n[1] * lab(ix,iy,iz-2).csf - n[0] * lab(ix,iy,iz-3).csf ),
                                  ifactor  * ( n[2] * lab(ix,iy,iz-1).csf - n[1] * lab(ix,iy,iz-2).csf ),
                                  ifactor  * (        lab(ix,iy,iz  ).csf - n[2] * lab(ix,iy,iz-1).csf ),
                                  ifactor  * ( n[3] * lab(ix,iy,iz+1).csf -        lab(ix,iy,iz  ).csf ),
                                  ifactor  * ( n[4] * lab(ix,iy,iz+2).csf - n[3] * lab(ix,iy,iz+1).csf ),
                                  ifactor  * ( n[5] * lab(ix,iy,iz+3).csf - n[4] * lab(ix,iy,iz+2).csf ),
                                  csf_flux_z );
                            
                            
#ifdef TissueMemory

                            /* In Z for tumour */
                            n[0] = lab(ix,iy,iz-3).p_w + lab(ix,iy,iz-3).p_g ;
                            n[1] = lab(ix,iy,iz-2).p_w + lab(ix,iy,iz-2).p_g ;
                            n[2] = lab(ix,iy,iz-1).p_w + lab(ix,iy,iz-1).p_g ;
                            n[3] = lab(ix,iy,iz+1).p_w + lab(ix,iy,iz+1).p_g ;
                            n[4] = lab(ix,iy,iz+2).p_w + lab(ix,iy,iz+2).p_g ;
                            n[5] = lab(ix,iy,iz+3).p_w + lab(ix,iy,iz+3).p_g ;
#else
                            n[0] = lab(ix,iy,iz-3).p_w + lab(ix,iy,iz-3).p_g + lab(ix,iy,iz-3).phi ;
                            n[1] = lab(ix,iy,iz-2).p_w + lab(ix,iy,iz-2).p_g + lab(ix,iy,iz-2).phi ;
                            n[2] = lab(ix,iy,iz-1).p_w + lab(ix,iy,iz-1).p_g + lab(ix,iy,iz-1).phi ;
                            n[3] = lab(ix,iy,iz+1).p_w + lab(ix,iy,iz+1).p_g + lab(ix,iy,iz+1).phi ;
                            n[4] = lab(ix,iy,iz+2).p_w + lab(ix,iy,iz+2).p_g + lab(ix,iy,iz+2).phi ;
                            n[5] = lab(ix,iy,iz+3).p_w + lab(ix,iy,iz+3).p_g + lab(ix,iy,iz+3).phi ;
                            
                            // normalise
                            n[0] = n[0] / (n[0] + lab(ix,iy,iz-3).p_csf );
                            n[1] = n[1] / (n[1] + lab(ix,iy,iz-2).p_csf );
                            n[2] = n[2] / (n[2] + lab(ix,iy,iz-1).p_csf );
                            n[3] = n[3] / (n[3] + lab(ix,iy,iz+1).p_csf );
                            n[4] = n[4] / (n[4] + lab(ix,iy,iz+2).p_csf );
                            n[5] = n[5] / (n[5] + lab(ix,iy,iz+3).p_csf );
#endif
                            
                            
#ifdef NoFluxBC
                            _applyNoFluxBC(n);
#endif
                            
                            // tumor
                            weno5( ifactor * ( n[1] * lab(ix,iy,iz-2).phi - n[0] * lab(ix,iy,iz-3).phi ),
                                  ifactor  * ( n[2] * lab(ix,iy,iz-1).phi - n[1] * lab(ix,iy,iz-2).phi ),
                                  ifactor  * (        lab(ix,iy,iz  ).phi - n[2] * lab(ix,iy,iz-1).phi ),
                                  ifactor  * ( n[3] * lab(ix,iy,iz+1).phi -        lab(ix,iy,iz  ).phi ),
                                  ifactor  * ( n[4] * lab(ix,iy,iz+2).phi - n[3] * lab(ix,iy,iz+1).phi ),
                                  ifactor  * ( n[5] * lab(ix,iy,iz+3).phi - n[4] * lab(ix,iy,iz+2).phi ),
                                  tumor_flux_z );
                            /*----*/
                            
                            // velocities
                            Real velocity_plus[3] = {
                                max(o(ix,iy,iz).ux, (Real) 0.),
                                max(o(ix,iy,iz).uy, (Real) 0.),
                                max(o(ix,iy,iz).uz, (Real) 0.)
                            };
                            
                            Real velocity_minus[3] = {
                                min(lab(ix,iy,iz).ux, (Real) 0.),
                                min(lab(ix,iy,iz).uy, (Real) 0.),
                                min(lab(ix,iy,iz).uz, (Real) 0.)
                                
                            };
                            
                            o(ix,iy,iz).dwmdt  = (-1.) * ( velocity_plus[0] * wm_flux_x[0] + velocity_minus[0] * wm_flux_x[1] +
                                                          velocity_plus[1] * wm_flux_y[0] + velocity_minus[1] * wm_flux_y[1] +
                                                          velocity_plus[2] * wm_flux_z[0] + velocity_minus[2] * wm_flux_z[1] );
                            
                            o(ix,iy,iz).dgmdt  = (-1.) * ( velocity_plus[0] * gm_flux_x[0] + velocity_minus[0] * gm_flux_x[1] +
                                                          velocity_plus[1] * gm_flux_y[0] + velocity_minus[1] * gm_flux_y[1] +
                                                          velocity_plus[2] * gm_flux_z[0] + velocity_minus[2] * gm_flux_z[1] );
                            
                            o(ix,iy,iz).dcsfdt = (-1.) * ( velocity_plus[0] * csf_flux_x[0] + velocity_minus[0] * csf_flux_x[1] +
                                                          velocity_plus[1] * csf_flux_y[0] + velocity_minus[1] * csf_flux_y[1] +
                                                          velocity_plus[2] * csf_flux_z[0] + velocity_minus[2] * csf_flux_z[1]);
                            
                            o(ix,iy,iz).dphidt = (-1.) * ( velocity_plus[0] * tumor_flux_x[0] + velocity_minus[0] * tumor_flux_x[1] +
                                                          velocity_plus[1] * tumor_flux_y[0] + velocity_minus[1] * tumor_flux_y[1] +
                                                          velocity_plus[2] * tumor_flux_z[0] + velocity_minus[2] * tumor_flux_z[1]);
                            
                        }
                        else
                        {
                            o(ix,iy,iz).dwmdt  = 0.;
                            o(ix,iy,iz).dgmdt  = 0.;
                            o(ix,iy,iz).dcsfdt = 0.;
                            o(ix,iy,iz).dphidt = 0.;
                            
                        }
                        
                    }
        }
    }
    
    
    inline void _applyNoFluxBC(Real (&n)[6] ) const
    {
        if( n[0] == 0. ){n[5] *= 2.; };
        if( n[1] == 0. ){n[4] *= 2.; };
        if( n[2] == 0. ){n[3] *= 2.; };
        if( n[3] == 0. ){n[2] *= 2.; };
        if( n[4] == 0. ){n[1] *= 2.; };
        if( n[5] == 0. ){n[0] *= 2.; };
    }
    
    
};


template<int nDim = 3>
struct TissueConvectivePart
{
    int stencil_start[3];
    int stencil_end[3];
    
    TissueConvectivePart()
    {
        stencil_start[0] = stencil_start[1] = -1;
        stencil_end[0]   = stencil_end[1]   = +2;
        stencil_start[2] = nDim==3 ? -1: 0;
        stencil_end[2]   = nDim==3 ? +2:+1;
    }
    
    TissueConvectivePart(const TissueConvectivePart&)
    {
        stencil_start[0] = stencil_start[1] = -1;
        stencil_end[0]   = stencil_end[1]   = +2;
        stencil_start[2] = nDim==3 ? -3: 0;
        stencil_end[2]   = nDim==3 ? +4:+1;
        
    }
    
    template<typename LabType, typename BlockType>
    inline void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        const double ih = 1./(info.h[0]);
        
        if(nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    if( lab(ix,iy).chi > 0)
                    {
                        Real uym = lab(ix  ,iy-1).uy;
                        Real uyp = lab(ix  ,iy+1).uy;
                        Real uy  = lab(ix,iy).uy;
                        
                        
                        Real uxm = lab(ix-1,iy  ).ux;
                        Real uxp = lab(ix+1,iy  ).ux;
                        Real ux = lab(ix,iy).ux;
                        
                        
                        _mean(uy,uym,uyp);
                        _mean(ux,uxm,uxp);
                        
                        
                        Real div_v = ih * (uxp - uxm + uyp - uym);
                        
                        o(ix,iy).dwmdt  -= lab(ix,iy).wm  * div_v;
                        o(ix,iy).dgmdt  -= lab(ix,iy).gm  * div_v;
                        o(ix,iy).dcsfdt -= lab(ix,iy).csf * div_v;
                        //o(ix,iy).dphidt -= lab(ix,iy).phi * div_v; // uncoment if tumor compression is considered
                    }
                    
                }
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        if( lab(ix,iy,iz).chi > 0)
                        {
                            Real uzm = lab(ix  ,iy  ,iz-1).uz;
                            Real uzp = lab(ix  ,iy  ,iz+1).uz;
                            Real uz  = lab(ix,iy,iz).uz;
                            
                            
                            Real uym = lab(ix  ,iy-1,iz  ).uy;
                            Real uyp = lab(ix  ,iy+1,iz  ).uy;
                            Real uy  = lab(ix,iy,iz).uy;
                            
                            
                            Real uxm = lab(ix-1,iy  ,iz  ).ux;
                            Real uxp = lab(ix+1,iy  ,iz  ).ux;
                            Real ux = lab(ix,iy,iz).ux;
                            
                            
                            _mean(uz,uzm,uzp);
                            _mean(uy,uym,uyp);
                            _mean(ux,uxm,uxp);
                            
                            
                            Real div_v = ih * (uxp - uxm + uyp - uym + uzp - uzm);
                            
                            o(ix,iy,iz).dwmdt  -= lab(ix,iy,iz).wm  * div_v;
                            o(ix,iy,iz).dgmdt  -= lab(ix,iy,iz).gm  * div_v;
                            o(ix,iy,iz).dcsfdt -= lab(ix,iy,iz).csf * div_v;
                            //o(ix,iy,iz).dphidt -= lab(ix,iy).phi * div_v;  // uncoment if tumor compression is considered
                        }
                        
                    }
        }
    }
    
    
    inline void _mean(Real u, Real& up, Real& um) const
    {
        up = 0.5 * (u + up);
        um = 0.5 * (u + um);
    }
    
};


template<int nDim = 3>
struct UpdateTissue
{
    double dt;
    
    UpdateTissue(double dt_)
    {
        dt = dt_;
    }
    
    UpdateTissue(const UpdateTissue&, double dt_ )
    {
        dt = dt_;
    }
    
    template<typename BlockType>
    inline void operator()(const BlockInfo& info, BlockType& o) const
    {
        if(nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    //update
                    o(ix, iy).wm    += dt * o(ix, iy).dwmdt;
                    o(ix, iy).gm    += dt * o(ix, iy).dgmdt;
                    o(ix, iy).csf   += dt * o(ix, iy).dcsfdt;
                    o(ix, iy).phi   += dt * o(ix, iy).dphidt;
                    
                    // reduce numerical error
                    o(ix,iy).wm     = max( o(ix,iy).wm,  (Real) 0.);
                    o(ix,iy).gm     = max( o(ix,iy).gm,  (Real) 0.);
                    o(ix,iy).csf    = max( o(ix,iy).csf, (Real) 0.);
                    o(ix,iy).phi    = max( o(ix,iy).phi, (Real) 0.);
                    o(ix,iy).phi    = min( o(ix,iy).phi, (Real) 1.);   // tumour is volume fraction, can't be more than 1
                    
                    
                    
                    /* recompute tissue percentage
                     If some tissue 100% tumor, do not modify it but keep what you remember about that tissue */
                    Real all = o(ix,iy).wm + o(ix,iy).gm + o(ix,iy).csf ;
                    
                    if( (all > 0. ) && (o(ix,iy).phi < 1.))
                    {
                        double i_all = 1./ all;
                        
                        o(ix,iy).p_w    = o(ix,iy).wm * i_all;
                        o(ix,iy).p_g    = o(ix,iy).gm * i_all;
                        o(ix,iy).p_csf  = o(ix,iy).csf * i_all;
                        
                        assert( o(ix,iy).p_w + o(ix,iy).p_g + o(ix,iy).p_csf < 1.1);
                    }
                }
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        //update
                        o(ix, iy, iz).wm    += dt * o(ix, iy, iz).dwmdt;
                        o(ix, iy, iz).gm    += dt * o(ix, iy, iz).dgmdt;
                        o(ix, iy, iz).csf   += dt * o(ix, iy, iz).dcsfdt;
                        o(ix, iy, iz).phi   += dt * o(ix, iy, iz).dphidt;
                        
                        // reduce numerical error
                        o(ix, iy, iz).wm     = max( o(ix, iy, iz).wm,  (Real) 0.);
                        o(ix, iy, iz).gm     = max( o(ix, iy, iz).gm,  (Real) 0.);
                        o(ix, iy, iz).csf    = max( o(ix, iy, iz).csf, (Real) 0.);
                        o(ix, iy, iz).phi    = max( o(ix, iy, iz).phi, (Real) 0.);
                        o(ix, iy, iz).phi    = min( o(ix, iy, iz).phi, (Real) 1.);   // tumour is volume fraction, can't be more than 1
                   
                        Real all = o(ix,iy,iz).wm + o(ix,iy,iz).gm + o(ix,iy,iz).csf ;

#ifndef TissueMemory
                        /* recompute tissue percentage
                         If some tissue 100% tumor, do not modify it but keep what you remember about that tissue */
                        if( (all > 0. ) && (o(ix,iy,iz).phi < 1.))
                        {
                            double i_all = 1./ all;
                            
                            o(ix,iy,iz).p_w    = o(ix,iy,iz).wm * i_all;
                            o(ix,iy,iz).p_g    = o(ix,iy,iz).gm * i_all;
                            o(ix,iy,iz).p_csf  = o(ix,iy,iz).csf * i_all;
                            
                            assert( o(ix,iy,iz).p_w + o(ix,iy,iz).p_g + o(ix,iy,iz).p_csf < 1.1);
                        }
#else
                        
    #ifdef CASE1
                        if( all > 0. )
                        {
                            double i_all = 1./ all;
                            
                            o(ix,iy,iz).p_w    = o(ix,iy,iz).wm * i_all;
                            o(ix,iy,iz).p_g    = o(ix,iy,iz).gm * i_all;
                            o(ix,iy,iz).p_csf  = o(ix,iy,iz).csf * i_all;
                            
                            assert( o(ix,iy,iz).p_w + o(ix,iy,iz).p_g + o(ix,iy,iz).p_csf < 1.1);
                        }
    #else
                        if( all > 0. )
                        {
                            double i_all = 1./ all;
                            
                            o(ix,iy,iz).p_w    = (1. - o(ix,iy,iz).phi ) * o(ix,iy,iz).wm * i_all;
                            o(ix,iy,iz).p_g    = (1. - o(ix,iy,iz).phi ) * o(ix,iy,iz).gm * i_all;
                            o(ix,iy,iz).p_csf  = (1. - o(ix,iy,iz).phi ) * o(ix,iy,iz).csf * i_all;
                            
                            
                            if(o(ix,iy,iz).p_w + o(ix,iy,iz).p_g + o(ix,iy,iz).p_csf + o(ix,iy,iz).phi > 1.1)
                            {
                                
                            printf("p_w=%f, p_g=%f, p_csf=%f, phi=%f \n", o(ix,iy,iz).p_w, o(ix,iy,iz).p_g, o(ix,iy,iz).p_csf, o(ix,iy,iz).phi );
                                abort();
                            }
                            //assert( o(ix,iy,iz).p_w + o(ix,iy,iz).p_g + o(ix,iy,iz).p_csf + o(ix,iy,iz).phi < 1.1);
                        }
    #endif
#endif
                        
                    }
        }
    }
};

#endif
