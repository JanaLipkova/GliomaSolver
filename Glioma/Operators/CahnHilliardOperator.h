//
//  CahnHilliardOperator.h
//  GliomaSolverXcode
//
//  Created by Lipkova on 21/09/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

/*Comput phase field function psi of binary domain defined by characteristic function

 Cahn-Hilliard equation:
 ----------------------------------------------------
 @ psi / @t = grad . ( M(psi) . grad(mu))  in Ω
 
 mu = g(psi)' - eps^2 * Δ psi      ; chemical potential
 g(psi)' = psi*(1-psi)*(0.5 - psi) ;double well potential
 M(psi) = psi*(1-psi) + tau        ; mobility

 + periodic BC			       in ∂Ω
 ------------------------------------------------
 tau = 0.001 (some small number, needed if psi(t=0) = {0,1})
 
 For details see 
  - section 2.2 in Tiegen et al; A DIFFUSE-INTERFACE APPROACH FOR MODELLING TRANSPORT, DIFFUSION AND ADSORPTION/DESORPTION OF MATERIAL QUANTITIES ON A DEFORMABLE INTERFACE
  - we assume zero velocity field
/*

 /* STEP 1: Compute Chemical Potential mu and define mobility term M(psi)  assuming double-well potential*/
template<int nDim = 3>
struct CahnHilliardPotential
{
    int stencil_start[3];
    int stencil_end[3];
    const int w;                 // width of smoothening
    
    CahnHilliardPotential(const int w_) :w(w_)
    {
        stencil_start[0] = stencil_start[1]= -1;
        stencil_end[0]   = stencil_end[1]  = +2;
        stencil_start[2] = nDim==3 ? -1:0;
        stencil_end[2]   = nDim==3?+2:+1;
    }
    
    CahnHilliardPotential(const CahnHilliardPotential& copy) :w(copy.w)
    {
        stencil_start[0] = stencil_start[1]= -1;
        stencil_end[0]   = stencil_end[1]  = +2;
        stencil_start[2] = nDim==3 ? -1:0;
        stencil_end[2]   = nDim==3?+2:+1;
    }
    
    template<typename LabType, typename BlockType>
    inline void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        double h		  = info.h[0];
        double ih2        = 1./(h*h);
        const double eps  = w*h;
        
        if (nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    const Real FluxInX = lab(ix-1,iy  ).pff + lab(ix+1,iy ).pff;
                    const Real FluxInY = lab(ix  ,iy-1).pff + lab(ix  ,iy+1).pff;
                    const Real FluxOut = 4.*lab(ix,iy).pff;
                    
                    const Real laplace = ih2 * (FluxInX + FluxInY - FluxOut);
                    
                    const Real tmp      = lab(ix,iy).pff * ( 1. - lab(ix,iy).pff) * (0.5 - lab(ix,iy).pff);  // derivative of potential g'(c)
                    const Real mobility = lab(ix,iy).pff * ( 1. - lab(ix,iy).pff);                   // mobility term M(c) = abs(c(1-c))
                    
                    o(ix,iy).mu  = tmp - eps*eps * laplace;
                    o(ix,iy).mob = (mobility < 0.) ? -mobility : mobility ;
                    
                }
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        const Real FluxInX = lab(ix-1,iy  ,iz  ).pff + lab(ix+1,iy  ,iz  ).pff;
                        const Real FluxInY = lab(ix  ,iy-1,iz  ).pff + lab(ix  ,iy+1,iz  ).pff;
                        const Real FluxInZ = lab(ix  ,iy  ,iz-1).pff + lab(ix  ,iy  ,iz+1).pff;
                        const Real FluxOut = 6.*lab(ix,iy,iz).pff;
                        
                        const Real laplace = ih2 * (FluxInX + FluxInY + FluxInZ - FluxOut);
                        
                        const Real tmp       = lab(ix,iy,iz).pff * ( 1. - lab(ix,iy,iz).pff) * (0.5 - lab(ix,iy,iz).pff);
                        const Real mobility  = lab(ix,iy,iz).pff * ( 1. - lab(ix,iy,iz).pff);
                        
                        o(ix,iy,iz).mu  = tmp - eps*eps * laplace;
                        o(ix,iy,iz).mob = (mobility < 0.) ? -mobility : mobility;
                    }
        }
        
    }
};



/* STEP 2: Compute phase field function */
template<int nDim = 3>
struct CahnHilliardEquation
{
    int stencil_start[3];
    int stencil_end[3];
    
    CahnHilliardEquation()
    {
        stencil_start[0] = stencil_start[1]= -1;
        stencil_end[0]   = stencil_end[1]  = +2;
        stencil_start[2] = nDim==3 ? -1:0;
        stencil_end[2]   = nDim==3?+2:+1;
    }
    
    CahnHilliardEquation(const CahnHilliardEquation& copy)
    {
        stencil_start[0] = stencil_start[1]= -1;
        stencil_end[0]   = stencil_end[1]  = +2;
        stencil_start[2] = nDim==3 ? -1:0;
        stencil_end[2]   = nDim==3?+2:+1;
    }
    
    template<typename LabType, typename BlockType>
    inline void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        double h		  = info.h[0];
        double ih2        = 1./(h*h);
        const double tau  = 1.e-7;
        
        if (nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    // compute mobility terms
                    Real mS = lab(ix  ,iy-1).mob + tau;
                    Real mN = lab(ix  ,iy+1).mob + tau;
                    Real mW = lab(ix-1,iy  ).mob + tau;
                    Real mE = lab(ix+1,iy  ).mob + tau;
                    Real m  = lab(ix  ,iy  ).mob + tau;
                    Real tmp1 = 0;
                    Real tmp2 = 0;
                    
                    // approximate intermidiet points
                    _mean(m, mS, mN, mW, mE, tmp1, tmp2);
                    //_harmonicAvg(m, mS, mN, mW, mE,tmp1,tmp2);
                    
                    Real fluxOutY = mS * lab(ix  ,iy-1).mu + mN * lab(ix  ,iy+1).mu;
                    Real fluxOutX = mW * lab(ix-1,iy  ).mu + mE * lab(ix+1,iy  ).mu;
                    Real fluxIn   = (mS + mN + mW + mE) * lab(ix,iy).mu;
                    
                    o(ix,iy).dpffdt = ih2 * (fluxOutX + fluxOutY - fluxIn );
                }
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        // compute mobility terms
                        Real mF =  lab(ix  ,iy  ,iz-1).mob + tau;
                        Real mB =  lab(ix  ,iy  ,iz+1).mob + tau;
                        Real mS =  lab(ix  ,iy-1,iz  ).mob + tau;
                        Real mN =  lab(ix  ,iy+1,iz  ).mob + tau;
                        Real mW =  lab(ix-1,iy  ,iz  ).mob + tau;
                        Real mE =  lab(ix+1,iy  ,iz  ).mob + tau;
                        Real m  =  lab(ix  ,iy  ,iz  ).mob + tau;
                        
                        // approximate intermidiet points
                        _mean(m, mS, mN, mW, mE, mB, mF);
                        
                        Real fluxOutZ = mF * lab(ix  ,iy  ,iz-1).mu + mB * lab(ix  ,iy  ,iz+1).mu;
                        Real fluxOutY = mS * lab(ix  ,iy-1,iz  ).mu + mN * lab(ix  ,iy+1,iz  ).mu;
                        Real fluxOutX = mW * lab(ix-1,iy  ,iz  ).mu + mE * lab(ix+1,iy  ,iz  ).mu;
                        Real fluxIn   = (mF + mB + mS + mN + mW + mE) * lab(ix,iy,iz).mu;
                        
                        o(ix,iy,iz).dpffdt = ih2 * (fluxOutZ + fluxOutY + fluxOutX - fluxIn );
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


/* STEP 3: Time integration with explicit Euler*/
template<int nDim = 3>
struct CahnHilliardUpdate
{
    double dt;
    
    CahnHilliardUpdate(double dt_) : dt(dt_)
    {    }
    
    CahnHilliardUpdate(const CahnHilliardUpdate& copy) : dt(copy.dt)
    {    }
    
    template<typename BlockType>
    inline void operator()(const BlockInfo& info, BlockType& o) const
    {
        if (nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    o(ix, iy).pff += dt * o(ix, iy).dpffdt;
                    o(ix, iy).pff  = (o(ix,iy).pff > 1.) ? 1. : o(ix,iy).pff;
                    o(ix, iy).pff  = (o(ix,iy).pff < 0.) ? 0. : o(ix,iy).pff;
                }
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        o(ix, iy, iz).pff += dt * o(ix, iy, iz).dpffdt;
                        o(ix, iy, iz).pff  = (o(ix,iy,iz).pff > 1.) ? 1. : o(ix,iy,iz).pff;
                        o(ix, iy, iz).pff  = (o(ix,iy,iz).pff < 0.) ? 0. : o(ix,iy,iz).pff;
                    }
        }
        
    }
};






