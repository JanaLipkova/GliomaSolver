//
//  Glioma_UQ_DataPreprocessing.h
//  GliomaSolverXcode
//
//  Created by Lipkova on 01/11/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

#ifndef __Glioma_UQ_DataPreprocessing__
#define __Glioma_UQ_DataPreprocessing__

#pragma once
#include "Glioma_Types.h"


class Glioma_UQ_DataPreprocessing: public Glioma
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
    bool                                    bVTK;
    string                                  PatientFileName;
    Real                                    cm[3];
    Real                                    volFLAIR;

    
    static void _ic(Grid<W,B>& grid, string PatientFileName);
    void        _normalisePET();
    void        _computeTumorStatistic();
    void        _dumpUQdata();
    void        _dumpInferenceROI();
    void		_dump(int counter);
    
public:
    Glioma_UQ_DataPreprocessing(int argc, const char ** argv);
    ~Glioma_UQ_DataPreprocessing();
    void run();
    
};


#endif /* defined(__Glioma_UQ_DataPreprocessing__) */
