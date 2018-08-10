//
//  PressureOperator.h
//  GliomaXcode
//
//  Created by Lipkova on 31/07/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

template<int nDim = 3>
struct PressureSourceOperator
{
    const Real rho;
    
    PressureSourceOperator(Real rho_):rho(rho_)
    { }
    
    PressureSourceOperator(const PressureSourceOperator& copy): rho(copy.rho)
    { }
    
    template<typename BlockType>
    inline void operator()(const BlockInfo& info, BlockType& o) const
    {
        if(nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                    o(ix, iy).f   = rho * o(ix, iy).phi * (1. - o(ix,iy).phi);
            
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                        o(ix, iy, iz).f   = rho * o(ix, iy,iz).phi * (1. - o(ix,iy,iz).phi) ;
            
        }
    }
};




template<int nDim = 3>
struct PressureGradient
{
    int stencil_start[3];
    int stencil_end[3];
    
    const bool bMobility;
    const Real Mwm, Mgm, Mcsf;    // mobility (i.e. hydraulic conductivity)
    
    PressureGradient(const bool bMobility_ = 0, const Real Mwm_ = 1, const Real Mgm_ = 1, const Real Mcsf_ = 1): bMobility(bMobility_), Mwm(Mwm_), Mgm(Mgm_), Mcsf(Mcsf_)
    {
        stencil_start[0]=stencil_start[1]= -1;
        stencil_end[0]  =stencil_end[1]  = +2;
        
        stencil_start[2] = nDim==3 ? -1:0;
        stencil_end[2]   = nDim==3 ? +2:+1;
    }
    
    PressureGradient(const PressureGradient& copy): bMobility(copy.bMobility), Mwm(copy.Mwm), Mgm(copy.Mgm), Mcsf(copy.Mcsf)
    {
        stencil_start[0]=stencil_start[1]= -1;
        stencil_end[0]  =stencil_end[1]  = +2;
        
        stencil_start[2] = nDim==3 ? -1:0;
        stencil_end[2]   = nDim==3 ? +2:+1;
    }
    
    
    template<typename LabType, typename BlockType>
    inline void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        
        const Real ih = 1./(info.h[0]);
        
        if(nDim==2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    if(lab(ix,iy).chi > 0)
                    {
                        Real pS = lab(ix,  iy-1).p;
                        Real pN = lab(ix,  iy+1).p;
                        Real pW = lab(ix-1,iy  ).p;
                        Real pE = lab(ix+1,iy  ).p;
                        Real p  = lab(ix  ,iy  ).p;
                        Real tmp1 = 0;
                        Real tmp2 = 0;
                        
                        // approximate intermidiet points
                        _mean(p, pW, pE, pS, pN, tmp1, tmp2);
                        
                        Real M = (bMobility) ? (Mwm * lab(ix,iy).p_w + Mgm * lab(ix,iy).p_g + Mcsf * lab(ix,iy).p_csf) : 1.;
                        
                        o(ix,iy).ux = -ih * M * (pE - pW) * lab(ix,iy).chi;
                        o(ix,iy).uy = -ih * M * (pN - pS) * lab(ix,iy).chi;
                    }
                    else
                    {
                        o(ix,iy).ux = 0.;
                        o(ix,iy).uy = 0.;
                    }
                }
            
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        if(lab(ix,iy,iz).chi > 0)
                        {
                            Real pB = lab(ix  ,iy  ,iz-1).p;
                            Real pF = lab(ix  ,iy  ,iz+1).p;
                            Real pS = lab(ix,  iy-1,iz  ).p;
                            Real pN = lab(ix,  iy+1,iz  ).p;
                            Real pW = lab(ix-1,iy  ,iz  ).p;
                            Real pE = lab(ix+1,iy  ,iz  ).p;
                            Real p  = lab(ix  ,iy  ,iz  ).p;
                            
                            // approximate intermidiet points
                            _mean(p, pW, pE, pS, pN, pB, pF);
                            
                            Real M = (bMobility) ? (Mwm * lab(ix,iy,iz).p_w + Mgm * lab(ix,iy,iz).p_g + Mcsf * lab(ix,iy,iz).p_csf) : 1.;
                            
                            o(ix,iy,iz).ux = -ih * M * (pE - pW) * lab(ix,iy,iz).chi;
                            o(ix,iy,iz).uy = -ih * M * (pN - pS) * lab(ix,iy,iz).chi;
                            o(ix,iy,iz).uz = -ih * M * (pF - pB) * lab(ix,iy,iz).chi;
                        }
                        else
                        {
                            o(ix,iy,iz).ux = 0.;
                            o(ix,iy,iz).uy = 0.;
                            o(ix,iy,iz).uz = 0.;
                        }
                    }
        }
        
    }
    
    
    inline void _mean(Real& m, Real& mS, Real& mN, Real& mW, Real& mE, Real& mB = 0., Real& mF = 0.) const
    {
        mS = 0.5 * (m + mS);
        mN = 0.5 * (m + mN);
        mW = 0.5 * (m + mW);
        mE = 0.5 * (m + mE);
        mB = 0.5 * (m + mB);
        mF = 0.5 * (m + mF);
    }
    
    
    inline void _harmonicAvg(Real& m, Real& mS, Real& mN, Real& mW, Real& mE, Real& mB = 0., Real& mF = 0.) const
    {
        mS = 2. * m * mS / (m + mS);
        mN = 2. * m * mN / (m + mN);
        mW = 2. * m * mW / (m + mW);
        mE = 2. * m * mE / (m + mE);
        mB = 2. * m * mB / (m + mB);
        mF = 2. * m * mF / (m + mF);
    }
};



/*
  Computes rhs of the pressure equation with time relaxation:
 //  ------------------------------------
 //    β ∂p / ∂t = ∇( M∇p) + ρ * φ(1-φ)
 //     + no flux BC
 //
 //     φ - tumor density
 //     ρ - tumor proliferation rate
 //     β - tissue compresibily
 //     M  - mobility (=hydraulic conductivity, i.e. easy with which pressure propaget through tissue)
 //  ------------------------------------
 */
template<int nDim = 3>
struct PressureTimeRelaxationOperator
{
    int stencil_start[3];
    int stencil_end[3];
    
    const Real mWM, mGM, mCSF, bWM, bGM, bCSF; // mobility and compressibility
    const Real rho; // tumor proliferation
    
    PressureTimeRelaxationOperator(const Real mCSF_, const Real mWM_, const Real mGM_, const Real bCSF_, const Real bWM_, const Real bGM_, const Real rho_): mCSF(mCSF_), mWM(mWM_), mGM(mGM_), bCSF(bCSF_), bWM(bWM_), bGM(bGM_), rho(rho_)
    {
        stencil_start[0] = stencil_start[1]= -1;
        stencil_end[0]   = stencil_end[1]  = +2;
        stencil_start[2] = nDim==3 ? -1: 0;
        stencil_end[2]   = nDim==3 ? +2:+1;
    }
    
    PressureTimeRelaxationOperator(const PressureTimeRelaxationOperator& copy): mCSF(copy.mCSF), mWM(copy.mWM), mGM(copy.mGM), bCSF(copy.bCSF), bWM(copy.bWM), bGM(copy.bGM), rho(copy.rho)
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
        
        Real mob[6];   // mobility coefficient
        Real chf[6];   // domain charact. func, chf=0 -> outside, chf=1 inside domain: use to apply BC
        Real mobLoc;   // to store mobility at current point (local)
        
        // K=buld modulus (defined as inverse of compressibility) to reduce amount od divisions
        const Real ibCSF = 1./bCSF;
        const Real ibWM =  1./bWM;
        const Real ibGM =  1./bGM;
        Real ib;
        
        if(nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    if( lab(ix,iy).chi > 0.)
                    {
                        ib     = lab(ix,iy).p_w * ibWM + lab(ix,iy).p_g * ibGM + lab(ix,iy).p_csf * ibCSF;
                        mobLoc = lab(ix,iy).p_w * mWM  + lab(ix,iy).p_g * mGM  + lab(ix,iy).p_csf * mCSF;
                      
                        mob[0]  = lab(ix-1,iy  ).p_w * mWM + lab(ix-1,iy  ).p_g * mGM + lab(ix-1,iy  ).p_csf * mCSF;
                        mob[1]  = lab(ix+1,iy  ).p_w * mWM + lab(ix+1,iy  ).p_g * mGM + lab(ix+1,iy  ).p_csf * mCSF;
                        mob[2]  = lab(ix  ,iy-1).p_w * mWM + lab(ix  ,iy-1).p_g * mGM + lab(ix  ,iy-1).p_csf * mCSF;
                        mob[3]  = lab(ix  ,iy+1).p_w * mWM + lab(ix  ,iy+1).p_g * mGM + lab(ix  ,iy+1).p_csf * mCSF;

                        _harmonic_mean(mob, mobLoc);
                        
                        chf[0] = lab(ix-1,iy  ).chi;
                        chf[1] = lab(ix+1,iy  ).chi;
                        chf[2] = lab(ix,  iy-1).chi;
                        chf[3] = lab(ix,  iy+1).chi;
                        
                        _applyNoFluxBC(mob, chf);
                        
                        // diffusion fluxes
                        double diffusionFluxIn  = ih2 * (mob[0]*lab(ix-1, iy).p +
                                                         mob[1]*lab(ix+1, iy).p +
                                                         mob[2]*lab(ix, iy-1).p +
                                                         mob[3]*lab(ix, iy+1).p  );
                        
                        double diffusionFluxOut = -( (mob[0] + mob[1] + mob[2] + mob[3]) * lab(ix, iy).p * ih2 );
                        double reactionFlux		= rho * lab(ix,iy).phi * ( 1. - lab(ix,iy).phi );
                        
                        o(ix, iy).dpdt =   ib * (diffusionFluxOut + diffusionFluxIn + reactionFlux) ;
                        
                    }
                    else
                        o(ix, iy).dpdt = 0.;
                    
                    
                }
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        if ( lab(ix,iy,iz).chi > 0.)
                        {
                            ib     = lab(ix,iy,iz).p_w * ibWM + lab(ix,iy,iz).p_g * ibGM + lab(ix,iy,iz).p_csf * ibCSF;
                            mobLoc = lab(ix,iy,iz).p_w * mWM  + lab(ix,iy,iz).p_g * mGM  + lab(ix,iy,iz).p_csf * mCSF;
                            
                            mob[0]  = lab(ix-1,iy  ,iz  ).p_w * mWM + lab(ix-1,iy  ,iz  ).p_g * mGM + lab(ix-1,iy  ,iz  ).p_csf * mCSF;
                            mob[1]  = lab(ix+1,iy  ,iz  ).p_w * mWM + lab(ix+1,iy  ,iz  ).p_g * mGM + lab(ix+1,iy  ,iz  ).p_csf * mCSF;
                            mob[2]  = lab(ix  ,iy-1,iz  ).p_w * mWM + lab(ix  ,iy-1,iz  ).p_g * mGM + lab(ix  ,iy-1,iz  ).p_csf * mCSF;
                            mob[3]  = lab(ix  ,iy+1,iz  ).p_w * mWM + lab(ix  ,iy+1,iz  ).p_g * mGM + lab(ix  ,iy+1,iz  ).p_csf * mCSF;
                            mob[4]  = lab(ix  ,iy  ,iz-1).p_w * mWM + lab(ix  ,iy  ,iz-1).p_g * mGM + lab(ix  ,iy  ,iz-1).p_csf * mCSF;
                            mob[5]  = lab(ix  ,iy  ,iz+1).p_w * mWM + lab(ix  ,iy  ,iz+1).p_g * mGM + lab(ix  ,iy  ,iz+1).p_csf * mCSF;
                            
                            _harmonic_mean(mob, mobLoc);
                            
                            chf[0] = lab(ix-1,iy  ,iz  ).chi;
                            chf[1] = lab(ix+1,iy  ,iz  ).chi;
                            chf[2] = lab(ix,  iy-1,iz  ).chi;
                            chf[3] = lab(ix,  iy+1,iz  ).chi;
                            chf[4] = lab(ix,  iy  ,iz-1).chi;
                            chf[5] = lab(ix,  iy  ,iz+1).chi;
                            
                            _applyNoFluxBC(mob, chf);
                            
                            // diffusion fluxes
                            double diffusionFluxIn  = ih2 * (mob[0]*lab(ix-1, iy, iz).p +
                                                             mob[1]*lab(ix+1, iy, iz).p +
                                                             mob[2]*lab(ix, iy-1, iz).p +
                                                             mob[3]*lab(ix, iy+1, iz).p +
                                                             mob[4]*lab(ix, iy, iz-1).p +
                                                             mob[5]*lab(ix, iy, iz+1).p   );
                            
                            double diffusionFluxOut = -( (mob[0] + mob[1] + mob[2] + mob[3] + mob[4] + mob[5]) * lab(ix, iy, iz).p * ih2 );
                            double reactionFlux		= rho * lab(ix,iy,iz).phi * ( 1. - lab(ix,iy,iz).phi );
                            
                            o(ix, iy, iz).dpdt =  ib * ( diffusionFluxOut + diffusionFluxIn + reactionFlux );
                            
                        }
                        else
                            o(ix, iy, iz).dpdt = 0.;
                        
                    }
            
        }
    }
    
    // Di,j = 2 * (Di * Dj / (Di + Dj)
    // set Di,j to zero if (Di + Dj = 0) i.e. no update and avoid division by zero
    inline void _harmonic_mean( Real (&m)[6], Real m_loc) const
    {
        Real eps = 1.0e-08; // to avoid divisin by zero
        
        m[0] = (m[0] + m_loc < eps) ? 0. : 2. * m[0] * m_loc / (m[0] + m_loc);
        m[1] = (m[1] + m_loc < eps) ? 0. : 2. * m[1] * m_loc / (m[1] + m_loc);
        m[2] = (m[2] + m_loc < eps) ? 0. : 2. * m[2] * m_loc / (m[2] + m_loc);
        m[3] = (m[3] + m_loc < eps) ? 0. : 2. * m[3] * m_loc / (m[3] + m_loc);
        
        if(nDim > 2)
        {
            m[4] = (m[4] + m_loc < eps) ? 0. : 2. * m[4] * m_loc / (m[4] + m_loc);
            m[5] = (m[5] + m_loc < eps) ? 0. : 2. * m[5] * m_loc / (m[5] + m_loc);
        }
        
    }
    
    inline void _applyNoFluxBC( Real (&mob)[6], Real chf[6] ) const
    {
        // n is domain char. func, use to apply bc by modifying the df term by the ghost point
        Real eps = 0.1;
        
        if(chf[0] < eps){mob[1] *= 2.0; }
        if(chf[1] < eps){mob[0] *= 2.0; }
        if(chf[2] < eps){mob[3] *= 2.0; }
        if(chf[3] < eps){mob[2] *= 2.0; }
        
        if(nDim > 2)
        {
            if(chf[4] < eps){mob[5] *= 2.0; }
            if(chf[5] < eps){mob[4] *= 2.0; }
        }
    }
    
};





