//
//  dat2VP.h
//  GliomaSolverXcode
//
//  Created by Lipkova on 31/08/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

#ifndef __GliomaSolverXcode__dat2VP__
#define __GliomaSolverXcode__dat2VP__

#pragma once
#include "Glioma_Types.h"
#include "../MRAG/MRAGio/DumpScalarToVP.h"


class dat2VP: public Glioma
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
    bool                                    isDone;
    bool                                    bVerbose;
    string                                  inputFileName;
    
    static void _ic(Grid<W,B>& grid, string inputFileName);
    void		_dumpVTK(int counter);
    void		_dumpVP(int counter);
    
    
public:
    dat2VP(int argc, const char ** argv);
    ~dat2VP();
    void run();
    
};

#endif /* defined(__GliomaSolverXcode__dat2VP__) */
