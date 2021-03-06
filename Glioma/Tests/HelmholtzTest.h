//
//  HelmholtzTest.h
//  GliomaSolver
//
//  Created by Lipkova on 31/03/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//


#pragma once
#include "../Glioma_Types.h"
#include "../Operators/HelmholtzSolver2D_Hypre.h"
#include "../Operators/HelmholtzSolver3D_Hypre.h"
#include "../Operators/HelmholtzSolver3D_Hypre_MPI.h"

class HelmholtzTest: public Glioma
{
private:
    Grid<W, B>								* grid;
    BlockProcessing							blockProcessing;
    Refiner_SpaceExtension					*refiner;
    Compressor								*compressor;
    //	BlockFWT<W, B, RD_Projector_Wavelets,true, 2>	blockfwt;  // use this if want refinment based on more channels
    BlockFWT<W, B, RD_Projector_Wavelets>	blockfwt;			// refinment based on single channel
    SpaceTimeSorter							stSorter;
    Profiler								profiler;
    ArgumentParser							parser;
    IO_VTK< W, B, RD_Projector_VTK >		vtk;
    BlockLab< B >							lab;
    int										numberOfIterations;
    bool                                    isDone;
    bool                                    bVerbose;
    bool                                    bProfiler;
    bool                                    bVTK;
    int                                     rank;
    int                                     nprocs;
    
    
    static void _ic_Square(Grid<W,B>& grid, int rank);
    void		_computeError();
    void        _dump(int counter);
        
public:
    HelmholtzTest(int argc, const char ** argv);
    ~HelmholtzTest();
    void run();
    
};
