//
//  Glioma_HG_UQ.h
//  GliomaXcode
//
//  Created by Lipkova on 10/06/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#ifndef _Glioma_ReactionDiffusion_
#define _Glioma_ReactionDiffusion_

#pragma once
#include "Glioma_Types.h"


class Glioma_ReactionDiffusion: public Glioma
{
private:
    Grid<W, B>								* grid;
    BlockProcessing							blockProcessing;
    Refiner_SpaceExtension					*refiner;
    Compressor								*compressor;
    BlockFWT<W, B, RD_Projector_Wavelets>	blockfwt;			// refinment based on single channel
    SpaceTimeSorter							stSorter;
    Profiler								profiler;
    ArgumentParser							parser;
    IO_VTK< W, B, RD_Projector_VTK >		vtk;
    BlockLab< B >							lab;
    int										numberOfIterations;
    double                                  whenToWrite;
    double                                  whenToWriteOffset;
    bool                                    isDone;
    bool                                    bAdaptivity;
    bool                                    bVerbose;
    bool                                    bVTK;
    bool                                    bUQ;
    string                                  PatientFileName;
    Real                                    L;
    Real                                    tumor_ic[3];
    
    static void _ic(Grid<W,B>& grid, string PatientFileName, Real& L, Real tumor_ic[3]);
    void        _reactionDiffusionStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, const Real Dw, const Real Dg, const Real rho,double dt);
    void		_dump(int counter);
    void        _dumpUQoutput();

public:
    Glioma_ReactionDiffusion(int argc, const char ** argv);
    ~Glioma_ReactionDiffusion();
    void run();
    
};

#endif /* defined(_Glioma_ReactionDiffusion_) */
