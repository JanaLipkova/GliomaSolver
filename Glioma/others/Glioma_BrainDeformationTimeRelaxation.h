//
//  Glioma_BrainDeformationTimeRelaxation.h
//  GliomaSolverXcode
//
//  Created by Lipkova on 08/08/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

#ifndef __Glioma_BrainDeformationTimeRelaxation__
#define __Glioma_BrainDeformationTimeRelaxation__

#pragma once
#include "Glioma_Types.h"


class Glioma_BrainDeformationTimeRelaxation: public Glioma
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
    bool                                    bProfiler;
    bool                                    bVTK;
    int                                     pID;
    Real                                    L;
    
    static void _ic(Grid<W,B>& grid, int pID, Real& L);
    static void _icSphere3Parts(Grid<W,B>& grid, Real& L );
    double      _estimate_dt(double dt, double CFL);
    Real        _compute_maxvel();
    void        _computeVelocities(BoundaryInfo* boundaryInfo ,const bool bMobility, std::vector<Real>* mobility);
    void        _pressureStep(BoundaryInfo* boundaryInfo, std::vector<Real> mobility, std::vector<Real> beta, const Real rho);
    void        _reactionDiffusionStep(BoundaryInfo* boundaryInfo, const Real Dw, const Real Dg, const Real rho);
    void        _advectionConvectionStep(BoundaryInfo* boundaryInfo);
    void        _timeUpdate(const int nParallelGranularity, const Real dt);
    void		_dump(int counter);
    
    
public:
    Glioma_BrainDeformationTimeRelaxation(int argc, const char ** argv);
    ~Glioma_BrainDeformationTimeRelaxation();
    void run();
    
};

#endif /* defined(__Glioma_BrainDeformationTimeRelaxation__) */
