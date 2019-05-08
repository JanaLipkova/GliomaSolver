//
//  TimeUpdateOperator
//  GliomaXcode
//
//  Created by Lipkova on 21/05/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#pragma once
#include "Glioma_Types.h"


template<int nDim = 3>
struct TimeUpdate
{
    double dt;
    
    TimeUpdate(double dt_)
    {
        dt = dt_;
    }
    
    TimeUpdate(const TimeUpdate&, double dt_ )
    {
        dt = dt_;
    }
    
    template<typename BlockType>
    inline void operator()(const BlockInfo& info, BlockType& o) const
    {
        Real eps = 1.e-06; // to avoid rounding errors
        
        if(nDim == 2)
        {
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    //time update
                    o(ix, iy).phi   += dt * o(ix, iy).dphidt;
                    o(ix, iy).wm    += dt * o(ix, iy).dwmdt;
                    o(ix, iy).gm    += dt * o(ix, iy).dgmdt;
                    o(ix, iy).csf   += dt * o(ix, iy).dcsfdt;
                    
                    // reduce numerical roundings errors
                    o(ix,iy).phi = ( o(ix,iy).phi > 1. ) ? 1. : o(ix,iy).phi;
                    o(ix,iy).phi = ( o(ix,iy).phi < 0. ) ? 0. : o(ix,iy).phi;
                    o(ix,iy).wm  = ( o(ix,iy).wm  < 0. ) ? 0. : o(ix,iy).wm;
                    o(ix,iy).gm  = ( o(ix,iy).gm  < 0. ) ? 0. : o(ix,iy).gm;
                    o(ix,iy).csf = ( o(ix,iy).csf < 0. ) ? 0. : o(ix,iy).csf;
                    
                    // compute tissue percantge
                    Real all = o(ix,iy).wm + o(ix,iy).gm + o(ix,iy).csf;
                    o(ix,iy).p_csf = (all > 0.) ? o(ix,iy).csf / all : 0.;
                    o(ix,iy).p_w   = (all > 0.) ? o(ix,iy).wm  / all : 0.;
                    o(ix,iy).p_g   = (all > 0.) ? o(ix,iy).gm  / all : 0.;
                }
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        //update
                        o(ix, iy, iz).phi   += dt * o(ix, iy, iz).dphidt;
                        o(ix, iy, iz).wm    += dt * o(ix, iy, iz).dwmdt;
                        o(ix, iy, iz).gm    += dt * o(ix, iy, iz).dgmdt;
                        o(ix, iy, iz).csf   += dt * o(ix, iy, iz).dcsfdt;
                        
                        // reduce numerical roundings errors
                        o(ix,iy,iz).phi = ( o(ix,iy,iz).phi > 1. ) ? 1. : o(ix,iy,iz).phi;
                        o(ix,iy,iz).phi = ( o(ix,iy,iz).phi < 0. ) ? 0. : o(ix,iy,iz).phi;
                        o(ix,iy,iz).wm  = ( o(ix,iy,iz).wm  < 0. ) ? 0. : o(ix,iy,iz).wm;
                        o(ix,iy,iz).gm  = ( o(ix,iy,iz).gm  < 0. ) ? 0. : o(ix,iy,iz).gm;
                        o(ix,iy,iz).csf = ( o(ix,iy,iz).csf < 0. ) ? 0. : o(ix,iy,iz).csf;
                        
                        // compute tissue percantge
                        Real all = o(ix,iy,iz).wm + o(ix,iy,iz).gm + o(ix,iy,iz).csf;
                        o(ix,iy,iz).p_csf = (all > 0.) ? o(ix,iy,iz).csf / all : 0.;
                        o(ix,iy,iz).p_w   = (all > 0.) ? o(ix,iy,iz).wm  / all : 0.;
                        o(ix,iy,iz).p_g   = (all > 0.) ? o(ix,iy,iz).gm  / all : 0.;
                        
                    }
        }
    }
};




template<int nDim = 3>
struct PressureTimeUpdate
{
    double dt;
    
    PressureTimeUpdate(double dt_)
    {
        dt = dt_;
    }
    
    PressureTimeUpdate(const PressureTimeUpdate&, double dt_ )
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
                    o(ix, iy).p += dt * o(ix, iy).dpdt;
                    o(ix, iy).p  = fabs(o(ix,iy).p);
                }
        }
        else
        {
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        o(ix, iy, iz).p += dt * o(ix, iy, iz).dpdt;
                        o(ix, iy, iz).p  = fabs(o(ix,iy, iz).p);
                    }
        }
    }
};
