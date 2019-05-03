//
//  ReactionDiffusionOperator.h
//  GliomaXcode
//
//  Created by Lipkova on 12/05/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//


/*
 Computes rhs of the reaction-diffusio equation for tumor:
 //  ------------------------------------
 //    ∂φ / ∂t = ∇( D∇φ) + ρ * φ(1-φ)
 //     + no flux BC
 //
 //     φ - tumor density
 //     ρ - tumor proliferation rate
 //     D - diffusivity
 //  ------------------------------------
 */
template<int nDim = 3>
struct ReactionDiffusionOperator
{
    int stencil_start[3];
    int stencil_end[3];
    
    const Real Dw, Dg, rho;
    
    ReactionDiffusionOperator(const Real Dw_, const Real Dg_, const Real rho_): Dw(Dw_), Dg(Dg_), rho(rho_)
    {
        stencil_start[0] = stencil_start[1]= -1;
        stencil_end[0]   = stencil_end[1]  = +2;
        stencil_start[2] = nDim==3 ? -1: 0;
        stencil_end[2]   = nDim==3 ? +2:+1;
    }
    
    ReactionDiffusionOperator(const ReactionDiffusionOperator& copy): Dw(copy.Dw), Dg(copy.Dg), rho(copy.rho)
    {
        stencil_start[0] = stencil_start[1]= -1;
        stencil_end[0]   = stencil_end[1]  = +2;
        stencil_start[2] = nDim==3 ? -1:0;
        stencil_end[2]   = nDim==3?+2:+1;
    }
    
    
    template<typename LabType, typename BlockType>
    inline void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        double h		= info.h[0];
        double ih2		= 1./(h*h);
        Real df[6];    // diffusion coefficient
        Real chf[6];   // domain charact. func, chf=0 -> outside, chf=1 inside domain: use to apply BC
        Real df_loc;   // diffusion at the current point (local)
        Real chf_loc;   // diffusion at the current point (local)

        if(nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    // check if we are in the brain domain = wm+gm + tumor
                    // need to include tumor for case of brain deformations, to deal with case where there is 100% tumor -> no healthy brain tissue but we are still in the brain
                    if ( lab(ix, iy).p_w + lab(ix, iy).p_g + lab(ix, iy).phi > 0.)
                    {
                        df_loc = lab(ix  ,iy  ).p_w * Dw + lab(ix  ,iy  ).p_g * Dg;
                        df[0]  = lab(ix-1,iy  ).p_w * Dw + lab(ix-1,iy  ).p_g * Dg;
                        df[1]  = lab(ix+1,iy  ).p_w * Dw + lab(ix+1,iy  ).p_g * Dg;
                        df[2]  = lab(ix  ,iy-1).p_w * Dw + lab(ix  ,iy-1).p_g * Dg;
                        df[3]  = lab(ix  ,iy+1).p_w * Dw + lab(ix  ,iy+1).p_g * Dg;
                        
                        _harmonic_mean(df, df_loc);

                        chf[0]  = lab(ix-1, iy  ).p_w + lab(ix-1, iy  ).p_w + lab(ix-1, iy  ).p_g ;
                        chf[1]  = lab(ix+1, iy  ).p_w + lab(ix+1, iy  ).p_w + lab(ix+1, iy  ).p_g ;
                        chf[2]  = lab(ix,   iy-1).p_w + lab(ix  , iy-1).p_w + lab(ix  , iy-1).p_g ;
                        chf[3]  = lab(ix,   iy+1).p_w + lab(ix  , iy+1).p_w + lab(ix  , iy+1).p_g ;
                        
                        _applyNoFluxBC(df,chf);
                        
                        // diffusion fluxes
                        double diffusionFluxIn  = ih2 * (df[0]*lab(ix-1, iy).phi +
                                                         df[1]*lab(ix+1, iy).phi +
                                                         df[2]*lab(ix, iy-1).phi +
                                                         df[3]*lab(ix, iy+1).phi  );
                        
                        double diffusionFluxOut = -( (df[0] + df[1] + df[2] + df[3]) * lab(ix, iy).phi * ih2 );
                        double reactionFlux		= rho * lab(ix,iy).phi * ( 1. - lab(ix,iy).phi );
                        
                        o(ix, iy).dphidt =   diffusionFluxOut + diffusionFluxIn + reactionFlux ;
                    }
                    else
                        o(ix, iy).dphidt = 0.;
                    
                }
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        // check if we are in the brain domain
                        if ( lab(ix,iy,iz).p_w + lab(ix,iy,iz).p_g + lab(ix,iy,iz).phi > 0.)
                        {
                            df_loc = lab(ix  ,iy  ,iz  ).p_w * Dw + lab(ix  ,iy  ,iz  ).p_g * Dg;
                            df[0]  = lab(ix-1,iy  ,iz  ).p_w * Dw + lab(ix-1,iy  ,iz  ).p_g * Dg;
                            df[1]  = lab(ix+1,iy  ,iz  ).p_w * Dw + lab(ix+1,iy  ,iz  ).p_g * Dg;
                            df[2]  = lab(ix  ,iy-1,iz  ).p_w * Dw + lab(ix  ,iy-1,iz  ).p_g * Dg;
                            df[3]  = lab(ix  ,iy+1,iz  ).p_w * Dw + lab(ix  ,iy+1,iz  ).p_g * Dg;
                            df[4]  = lab(ix  ,iy  ,iz-1).p_w * Dw + lab(ix  ,iy  ,iz-1).p_g * Dg;
                            df[5]  = lab(ix  ,iy  ,iz+1).p_w * Dw + lab(ix  ,iy  ,iz+1).p_g * Dg;
                    
                           _harmonic_mean(df, df_loc);

                            chf[0] = lab(ix-1,iy,  iz  ).phi + lab(ix-1,iy,  iz  ).p_w + lab(ix-1,iy,  iz  ).p_g ;
                            chf[1] = lab(ix+1,iy,  iz  ).phi + lab(ix+1,iy,  iz  ).p_w + lab(ix+1,iy,  iz  ).p_g ;
                            chf[2] = lab(ix  ,iy-1,iz  ).phi + lab(ix  ,iy-1,iz  ).p_w + lab(ix  ,iy-1,iz  ).p_g ;
                            chf[3] = lab(ix  ,iy+1,iz  ).phi + lab(ix  ,iy+1,iz  ).p_w + lab(ix  ,iy+1,iz  ).p_g ;
                            chf[4] = lab(ix  ,iy,  iz-1).phi + lab(ix  ,iy,  iz-1).p_w + lab(ix  ,iy,  iz-1).p_g ;
                            chf[5] = lab(ix  ,iy,  iz+1).phi + lab(ix  ,iy,  iz+1).p_w + lab(ix  ,iy,  iz+1).p_g ;
                            
                            _applyNoFluxBC(df,chf);
                            
                            // diffusion fluxes
                            double diffusionFluxIn  = ih2 * (df[0]*lab(ix-1, iy, iz).phi +
                                                             df[1]*lab(ix+1, iy, iz).phi +
                                                             df[2]*lab(ix, iy-1, iz).phi +
                                                             df[3]*lab(ix, iy+1, iz).phi +
                                                             df[4]*lab(ix, iy, iz-1).phi +
                                                             df[5]*lab(ix, iy, iz+1).phi   );
                            
                            double diffusionFluxOut = -( (df[0] + df[1] + df[2] + df[3] + df[4] + df[5]) * lab(ix, iy, iz).phi * ih2 );
                            double reactionFlux		= rho * lab(ix,iy,iz).phi * ( 1. - lab(ix,iy,iz).phi );
                            
                            o(ix, iy, iz).dphidt =  diffusionFluxOut + diffusionFluxIn + reactionFlux ;
                        }
                        else
                            o(ix, iy, iz).dphidt = 0.;
                    }
            
        }
    }
    
    
    
    // Di,j = 2 * (Di * Dj / (Di + Dj)
    // set Di,j to zero if (Di + Dj = 0) i.e. no update and avoid division by zero
    inline void _harmonic_mean( Real (&df)[6], Real df_loc) const
    {
        Real eps = 1.0e-08; // to avoid divisin by zero
        
        df[0] = (df[0] + df_loc < eps) ? 0. : 2. * df[0] * df_loc / (df[0] + df_loc );
        df[1] = (df[1] + df_loc < eps) ? 0. : 2. * df[1] * df_loc / (df[1] + df_loc );
        df[2] = (df[2] + df_loc < eps) ? 0. : 2. * df[2] * df_loc / (df[2] + df_loc );
        df[3] = (df[3] + df_loc < eps) ? 0. : 2. * df[3] * df_loc / (df[3] + df_loc );
        
        if(nDim > 2)
        {
            df[4] = (df[4] + df_loc < eps) ? 0. : 2. * df[4] * df_loc / (df[4]  + df_loc );
            df[5] = (df[5] + df_loc < eps) ? 0. : 2. * df[5] * df_loc / (df[5]  + df_loc );
        }
        
    }
    
    
    inline void _applyNoFluxBC( Real (&df)[6], Real n[6] ) const
    {
        // n is domain char. func, use to apply bc by modifying the df term by the ghost point
        Real eps = 0.1;
        
        if(n[0] < eps){df[1] *= 2.0; }
        if(n[1] < eps){df[0] *= 2.0; }
        if(n[2] < eps){df[3] *= 2.0; }
        if(n[3] < eps){df[2] *= 2.0; }
        
        if(nDim > 2)
        {
            if(n[4] < eps){df[5] *= 2.0; }
            if(n[5] < eps){df[4] *= 2.0; }
        }
    }
    
};




template<int nDim = 3>
struct UpdateTumor
{
    double dt;
    
    UpdateTumor(double dt_):dt(dt_)
    { }
    
    UpdateTumor(const UpdateTumor& copy):dt(copy.dt)
    { }
    
    template<typename BlockType>
    inline void operator()(const BlockInfo& info, BlockType& o) const
    {
        if(nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    o(ix, iy).phi += dt * o(ix, iy).dphidt;
                    o(ix, iy).phi = max((Real)0., o(ix,iy).phi);
                    o(ix, iy).phi = min((Real)1., o(ix, iy).phi);
                }
            
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        o(ix, iy, iz).phi += dt * o(ix, iy, iz).dphidt;
                        o(ix, iy, iz).phi = max((Real)0., o(ix,iy,iz).phi);
                        o(ix, iy, iz).phi = min((Real)1., o(ix,iy,iz).phi);
                    }
            
        }
        
    }
};

