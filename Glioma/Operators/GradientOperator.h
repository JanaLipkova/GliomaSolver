/*
 *  GradientOperator.h
 *  GliomaXcode
 *
 *  Created by Lipkova on 10/6/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */


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






