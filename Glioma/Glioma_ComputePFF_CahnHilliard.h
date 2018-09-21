//
//  Glioma_ComputePFF_CahnHilliard.h
//  GliomaSolverXcode
//
//  Created by Lipkova on 21/09/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

#pragma once
#include "Glioma_Types.h"

class Glioma_ComputePFF_CahnHilliard: public Glioma
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
    bool                                    bAdaptivity;
    bool                                    bProfiler;
    bool                                    bVerbose;
    bool                                    bVTK;
    string                                  PatientFileName;
    Real                                    L;
    bool                                    isDone;

    static void _ic(Grid<W,B>& grid, string PatientFileName, Real& L);
    
    void        _CahnHilliardStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, double dt, const int w);
    void		_dump(int counter);
    void        _dumpBinary(int counter);
    
public:
    Glioma_ComputePFF_CahnHilliard(int argc, const char ** argv);
    ~Glioma_ComputePFF_CahnHilliard();
    void run();
    
};