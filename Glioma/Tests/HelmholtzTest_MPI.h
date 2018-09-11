//
//  HelmholtzTest_MPI.h
//  GliomaSolver
//
//  Created by Lipkova on 11/09/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

#pragma once
#include "../Glioma_Types.h"
#include "../Operators/HelmholtzSolver3D_Hypre_MPI.h"

class HelmholtzTest_MPI: public Glioma
{
private:
    Grid<W, B>								* grid;
    BlockProcessing							blockProcessing;
    Refiner_SpaceExtension					*refiner;
    Compressor								*compressor;
    BlockFWT<W, B, RD_Projector_Wavelets>	blockfwt;			
    SpaceTimeSorter							stSorter;
    Profiler								profiler;
    ArgumentParser							parser;
    IO_VTK< W, B, RD_Projector_VTK >		vtk;
    BlockLab< B >							lab;
    bool                                    isDone;
    bool                                    bVerbose;
    bool                                    bProfiler;
    bool                                    bVTK;
    
    static void _ic_Square(Grid<W,B>& grid);
    void		_computeError();
    void        _dump(int counter);
    
    HelmholtzSolver3D_Hypre_MPI helmholtz_solver3D_MPI;
    
public:
    HelmholtzTest_MPI(int argc, const char ** argv);
    ~HelmholtzTest_MPI();
    void run();
    
};
