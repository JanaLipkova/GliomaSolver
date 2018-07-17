//
//  Glioma_ReactionDiffusionOperator.h
//  GliomaXcode
//
//  Created by Lipkova on 12/05/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#ifndef GliomaXcode_Glioma_ReactionDiffusionOperator_h
#define GliomaXcode_Glioma_ReactionDiffusionOperator_h

template<int nDim = 3>
struct Glioma_ReactionDiffusionOperator
{
    int stencil_start[3];
    int stencil_end[3];
    
    const Real Dw, Dg, rho;
    
    Glioma_ReactionDiffusionOperator(const Real Dw_, const Real Dg_, const Real rho_): Dw(Dw_), Dg(Dg_), rho(rho_)
    {
        stencil_start[0] = stencil_start[1]= -1;
        stencil_end[0]   = stencil_end[1]  = +2;
        stencil_start[2] = nDim==3 ? -1: 0;
        stencil_end[2]   = nDim==3 ? +2:+1;
    }
    
    Glioma_ReactionDiffusionOperator(const Glioma_ReactionDiffusionOperator& copy): Dw(copy.Dw), Dg(copy.Dg), rho(copy.rho)
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
        double df[6];
        
        
        
        if(nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    df[0] = 0.0; df[1] = 0.0; df[2] = 0.0; df[3] = 0.0;
                    
                    // check if we are in the brain domain
                    if ( (lab(ix, iy).p_w > 0.0) || (lab(ix, iy).p_g > 0.0) )
                    {
                        
                        // Harmonic averages of piecewise constant diffusion coefficients
                        // Bernstein 2005 is wrong: factor of 2 is missing
                        // TRICK: double x = 1./0; then double y = 1./x is 0 what is cooool
                        // so we directly obtain diffusion zero for gp out of domain :)))
                        
                        df[0] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy).p_w*Dw + lab(ix, iy).p_g*Dg) ) + (1.0 / (lab(ix-1, iy).p_w*Dw + lab(ix-1, iy).p_g*Dg) ) ) );
                        df[1] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy).p_w*Dw + lab(ix, iy).p_g*Dg) ) + (1.0 / (lab(ix+1, iy).p_w*Dw + lab(ix+1, iy).p_g*Dg) ) ) );
                        df[2] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy).p_w*Dw + lab(ix, iy).p_g*Dg) ) + (1.0 / (lab(ix, iy-1).p_w*Dw + lab(ix, iy-1).p_g*Dg) ) ) );
                        df[3] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy).p_w*Dw + lab(ix, iy).p_g*Dg) ) + (1.0 / (lab(ix, iy+1).p_w*Dw + lab(ix, iy+1).p_g*Dg) ) ) );
                        
                        
                    }
                    
                    // Neumann no flux boundary condition, 2nd order using ghosts
                    // need correction by factor of 2
                    // if some df = 0, means it is boundary point,the oposit direction need to be multiply by 2
                    if ( (df[0] == 0) && ( lab(ix-1,iy).phi == 0) ){ df[1] *= 2.0; }
                    if ( (df[1] == 0) && ( lab(ix+1,iy).phi == 0) ){ df[0] *= 2.0; }
                    if ( (df[2] == 0) && ( lab(ix,iy-1).phi == 0) ){ df[3] *= 2.0; }
                    if ( (df[3] == 0) && ( lab(ix,iy+1).phi == 0) ){ df[2] *= 2.0; }
                    
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
            
#define Case2
            
#ifdef TissueMemory
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        df[0] = 0.0; df[1] = 0.0; df[2] = 0.0; df[3] = 0.0; df[4] = 0.0; df[5] = 0.0;
                        
                        // check if we are in the brain domain
                        if ( (lab(ix, iy, iz).p_w > 0.0) || (lab(ix, iy, iz).p_g > 0.0) )
                        {
                            // Harmonic averages of piecewise constant diffusion coefficients
                            // Bernstein 2005 is wrong: factor of 2 is missing
                            // TRICK: double x = 1./0; then double y = 1./x is 0 what is cooool
                            // so we directly obtain diffusion zero for gp out of domain :)))
                            
                            
                            df[0] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix-1, iy, iz).p_w*Dw + lab(ix-1, iy, iz).p_g*Dg) ) ) );
                            df[1] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix+1, iy, iz).p_w*Dw + lab(ix+1, iy, iz).p_g*Dg) ) ) );
                            df[2] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix, iy-1, iz).p_w*Dw + lab(ix, iy-1, iz).p_g*Dg) ) ) );
                            df[3] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix, iy+1, iz).p_w*Dw + lab(ix, iy+1, iz).p_g*Dg) ) ) );
                            df[4] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix, iy, iz-1).p_w*Dw + lab(ix, iy, iz-1).p_g*Dg) ) ) );
                            df[5] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix, iy, iz+1).p_w*Dw + lab(ix, iy, iz+1).p_g*Dg) ) ) );
                            
                        }
                        
                        // Neumann no flux boundary condition, 2nd order using ghosts
                        // need correction by factor of 2
                        // if some df = 0, means it is boundary point,the oposit direction need to be multiply by 2
                        if ( df[0] == 0 ){ df[1] *= 2.0; }
                        if ( df[1] == 0 ){ df[0] *= 2.0; }
                        if ( df[2] == 0 ){ df[3] *= 2.0; }
                        if ( df[3] == 0 ){ df[2] *= 2.0; }
                        if ( df[4] == 0 ){ df[5] *= 2.0; }
                        if ( df[5] == 0 ){ df[4] *= 2.0; }
                        
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
            
#else
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
                            // Bernstein 2005 is wrong: factor of 2 is missing
                            // TRICK: double x = 1./0; then double y = 1./x is 0 what is cooool
                            // so we directly obtain diffusion zero for gp out of domain :)))
                            
#ifdef Case1
                            df[0] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix-1, iy, iz).p_w*Dw + lab(ix-1, iy, iz).p_g*Dg) ) ) );
                            df[1] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix+1, iy, iz).p_w*Dw + lab(ix+1, iy, iz).p_g*Dg) ) ) );
                            df[2] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix, iy-1, iz).p_w*Dw + lab(ix, iy-1, iz).p_g*Dg) ) ) );
                            df[3] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix, iy+1, iz).p_w*Dw + lab(ix, iy+1, iz).p_g*Dg) ) ) );
                            df[4] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix, iy, iz-1).p_w*Dw + lab(ix, iy, iz-1).p_g*Dg) ) ) );
                            df[5] = 2.0*(1.0 / ( (1.0 / (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) ) + (1.0 / (lab(ix, iy, iz+1).p_w*Dw + lab(ix, iy, iz+1).p_g*Dg) ) ) );
#endif
                            
#ifdef Case2
                            Real brainW = lab(ix-1,iy,  iz  ).p_w + lab(ix-1,iy,  iz  ).p_g;
                            Real brain0 = lab(ix,  iy,  iz  ).p_w + lab(ix,  iy,  iz  ).p_g;
                            Real brainE = lab(ix+1,iy,  iz  ).p_w + lab(ix+1,iy,  iz  ).p_g;
                            Real brainS = lab(ix,  iy-1,iz  ).p_w + lab(ix,  iy-1,iz  ).p_g;
                            Real brainN = lab(ix,  iy+1,iz  ).p_w + lab(ix,  iy+1,iz  ).p_g;
                            Real brainB = lab(ix,  iy,  iz-1).p_w + lab(ix,  iy,  iz-1).p_g;
                            Real brainF = lab(ix,  iy,  iz+1).p_w + lab(ix,  iy,  iz+1).p_g;
                            
                            
                            df[0] = 2.0*(1.0 / ( (1.0 / ( (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) / brain0) ) + (1.0 / ( (lab(ix-1, iy, iz).p_w*Dw + lab(ix-1, iy, iz).p_g*Dg) / brainW ) ) ) );
                            df[1] = 2.0*(1.0 / ( (1.0 / ( (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) / brain0) ) + (1.0 / ( (lab(ix+1, iy, iz).p_w*Dw + lab(ix+1, iy, iz).p_g*Dg) / brainE ) ) ) );
                            df[2] = 2.0*(1.0 / ( (1.0 / ( (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) / brain0) ) + (1.0 / ( (lab(ix, iy-1, iz).p_w*Dw + lab(ix, iy-1, iz).p_g*Dg) / brainS ) ) ) );
                            df[3] = 2.0*(1.0 / ( (1.0 / ( (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) / brain0) ) + (1.0 / ( (lab(ix, iy+1, iz).p_w*Dw + lab(ix, iy+1, iz).p_g*Dg) / brainN ) ) ) );
                            df[4] = 2.0*(1.0 / ( (1.0 / ( (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) / brain0) ) + (1.0 / ( (lab(ix, iy, iz-1).p_w*Dw + lab(ix, iy, iz-1).p_g*Dg) / brainB ) ) ) );
                            df[5] = 2.0*(1.0 / ( (1.0 / ( (lab(ix, iy, iz).p_w*Dw + lab(ix, iy, iz).p_g*Dg) / brain0) ) + (1.0 / ( (lab(ix, iy, iz+1).p_w*Dw + lab(ix, iy, iz+1).p_g*Dg) / brainF ) ) ) );
#endif
                            
                            
                        }
                        
                        // Neumann no flux boundary condition, 2nd order using ghosts
                        // need correction by factor of 2
                        // if some df = 0, means it is boundary point,the oposit direction need to be multiply by 2
                        // Apply BC only if there is no tumour
                        if ( df[0] == 0 && lab(ix-1,iy,  iz  ).phi > 0. ){ df[1] *= 2.0; }
                        if ( df[1] == 0 && lab(ix+1,iy,  iz  ).phi > 0. ){ df[0] *= 2.0; }
                        if ( df[2] == 0 && lab(ix,  iy-1,iz  ).phi > 0. ){ df[3] *= 2.0; }
                        if ( df[3] == 0 && lab(ix,  iy+1,iz  ).phi > 0. ){ df[2] *= 2.0; }
                        if ( df[4] == 0 && lab(ix,  iy,  iz-1).phi > 0. ){ df[5] *= 2.0; }
                        if ( df[5] == 0 && lab(ix,  iy,  iz+1).phi > 0. ){ df[4] *= 2.0; }
                        
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
#endif
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



template<int nDim = 3>
struct UpdatePressureSource
{
    const Real rho;
    
    UpdatePressureSource(Real rho_):rho(rho_)
    { }
    
    UpdatePressureSource(const UpdatePressureSource& copy): rho(copy.rho)
    { }
    
    template<typename BlockType>
    inline void operator()(const BlockInfo& info, BlockType& o) const
    {
        double tau = 1.e-10; // phase field f. correction, pff = max(pff,tau)
        
        if(nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    Real pff      = ( o(ix,iy).psi > tau ) ? o(ix,iy).psi : 0.  ;
                    o(ix, iy).f   = rho * o(ix, iy).phi * (1. - o(ix,iy).phi) * pff ;
                }
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        Real pff          = ( o(ix,iy,iz).psi > tau ) ? o(ix,iy,iz).psi : 0.  ;
                        o(ix, iy, iz).f   = rho * o(ix, iy,iz).phi * (1. - o(ix,iy,iz).phi) * pff ;
                    }
        }
        
    }
};

#endif


