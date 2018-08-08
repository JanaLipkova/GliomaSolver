//
//  ReactionDiffusionOperator.h
//  GliomaXcode
//
//  Created by Lipkova on 12/05/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//


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
        Real df[6];  // diffusion coefficient
        Real n[6];   // domain charact. func, n=0 -> outside, n=1 inside domain: use to apply BC
        
        if(nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    df[0] = 0.0; df[1] = 0.0; df[2] = 0.0; df[3] = 0.0;
                    
                    // check if we are in the brain domain = wm+gm + tumor
                    // need to include tumor for case of brain deformations, to deal with case where there is 100% tumor -> no healthy brain tissue but we are still in the brain
                    Real tissue = lab(ix, iy).p_w + lab(ix, iy).p_g + lab(ix, iy).phi;
                    if ( tissue > 0.)
                    {
                        // Harmonic averages of piecewise constant diffusion coefficients
                        // TRICK: double x = 1./0; then double y = 1./x is 0 what is cooool
                        // so we directly obtain diffusion zero for gp out of domain :)))
                        df[0] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy).p_w*Dw + lab(ix, iy).p_g*Dg) ) + (1.0 / (lab(ix-1, iy).p_w*Dw + lab(ix-1, iy).p_g*Dg) ) ) );
                        df[1] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy).p_w*Dw + lab(ix, iy).p_g*Dg) ) + (1.0 / (lab(ix+1, iy).p_w*Dw + lab(ix+1, iy).p_g*Dg) ) ) );
                        df[2] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy).p_w*Dw + lab(ix, iy).p_g*Dg) ) + (1.0 / (lab(ix, iy-1).p_w*Dw + lab(ix, iy-1).p_g*Dg) ) ) );
                        df[3] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy).p_w*Dw + lab(ix, iy).p_g*Dg) ) + (1.0 / (lab(ix, iy+1).p_w*Dw + lab(ix, iy+1).p_g*Dg) ) ) );
                    }
                    
                    // include phi to n to account for the case where u=1, but tissue=0 due to deformation
                    n[0] = lab(ix-1, iy  ).p_w + lab(ix-1, iy  ).p_g + lab(ix-1, iy  ).phi;
                    n[1] = lab(ix+1, iy  ).p_w + lab(ix+1, iy  ).p_g + lab(ix+1, iy  ).phi;
                    n[2] = lab(ix  , iy-1).p_w + lab(ix  , iy-1).p_g + lab(ix  , iy-1).phi;
                    n[3] = lab(ix  , iy+1).p_w + lab(ix  , iy+1).p_g + lab(ix  , iy+1).phi;
                    
                    _applyNoFluxBC(df,n);
                    
                    // diffusion fluxes
                    double diffusionFluxIn  = ih2 * (df[0]*lab(ix-1, iy).phi +
                                                     df[1]*lab(ix+1, iy).phi +
                                                     df[2]*lab(ix, iy-1).phi +
                                                     df[3]*lab(ix, iy+1).phi  );
                    
                    double diffusionFluxOut = -( (df[0] + df[1] + df[2] + df[3]) * lab(ix, iy).phi * ih2 );
                    double reactionFlux		= rho * lab(ix,iy).phi * ( 1. - lab(ix,iy).phi );
                    
                    o(ix, iy).dphidt =   diffusionFluxOut + diffusionFluxIn + reactionFlux ;
                    
                }
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        df[0] = 0.0; df[1] = 0.0; df[2] = 0.0; df[3] = 0.0; df[4] = 0.0; df[5] = 0.0;
                        
                        // check if we are in the brain domain
                        Real tissue = lab(ix,iy,iz).p_w + lab(ix,iy,iz).p_g + lab(ix,iy,iz).phi;
                        if ( tissue > 0.)
                        {
                            // Harmonic averages of piecewise constant diffusion coefficients
                            // TRICK: double x = 1./0; then double y = 1./x is 0 what is cooool
                            // so we directly obtain diffusion zero for gp out of domain :)))
                            df[0] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix-1, iy, iz).p_w*Dw + lab(ix-1, iy, iz).p_g*Dg) ) ) );
                            df[1] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix+1, iy, iz).p_w*Dw + lab(ix+1, iy, iz).p_g*Dg) ) ) );
                            df[2] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix, iy-1, iz).p_w*Dw + lab(ix, iy-1, iz).p_g*Dg) ) ) );
                            df[3] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix, iy+1, iz).p_w*Dw + lab(ix, iy+1, iz).p_g*Dg) ) ) );
                            df[4] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix, iy, iz-1).p_w*Dw + lab(ix, iy, iz-1).p_g*Dg) ) ) );
                            df[5] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix, iy, iz+1).p_w*Dw + lab(ix, iy, iz+1).p_g*Dg) ) ) );

                        }
                        
                        n[0] = lab(ix-1,iy,  iz  ).p_w + lab(ix-1,iy,  iz  ).p_g + lab(ix-1,iy,  iz  ).phi;
                        n[1] = lab(ix+1,iy,  iz  ).p_w + lab(ix+1,iy,  iz  ).p_g + lab(ix+1,iy,  iz  ).phi;
                        n[2] = lab(ix  ,iy-1,iz  ).p_w + lab(ix  ,iy-1,iz  ).p_g + lab(ix  ,iy-1,iz  ).phi;
                        n[3] = lab(ix  ,iy+1,iz  ).p_w + lab(ix  ,iy+1,iz  ).p_g + lab(ix  ,iy+1,iz  ).phi;
                        n[4] = lab(ix  ,iy,  iz-1).p_w + lab(ix  ,iy,  iz-1).p_g + lab(ix  ,iy,  iz-1).phi;
                        n[5] = lab(ix  ,iy,  iz+1).p_w + lab(ix  ,iy,  iz+1).p_g + lab(ix  ,iy,  iz+1).phi;

                        _applyNoFluxBC(df,n);
                        
                        // diffusion fluxes
                        double diffusionFluxIn  = ih2 * (df[0]*lab(ix-1, iy, iz).phi +
                                                         df[1]*lab(ix+1, iy, iz).phi +
                                                         df[2]*lab(ix, iy-1, iz).phi +
                                                         df[3]*lab(ix, iy+1, iz).phi +
                                                         df[4]*lab(ix, iy, iz-1).phi +
                                                         df[5]*lab(ix, iy, iz+1).phi   );
                        
                        double diffusionFluxOut = -( (df[0] + df[1] + df[2] + df[3] + df[4] + df[5]) * lab(ix, iy, iz).phi * ih2 );
                        double reactionFlux		= rho * lab(ix,iy,iz).phi * ( 1. - lab(ix,iy,iz).phi );
                        
                        o(ix, iy, iz).dphidt =   diffusionFluxOut + diffusionFluxIn + reactionFlux ;
                        
                    }

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

