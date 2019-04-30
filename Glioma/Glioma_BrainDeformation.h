//
//  Glioma_BrainDeformation.h
//  GliomaSolver
//
//  Created by Lipkova on 17/07/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

#ifndef __Glioma_BrainDeformation__
#define __Glioma_BrainDeformation__


#pragma once
#include "Glioma_Types.h"
#include "../Operators/HelmholtzSolver3D_Hypre.h"
#include "Operators/HelmholtzSolver3D_Hypre_MPI.h"

class Glioma_BrainDeformation: public Glioma
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
    bool                                    bAllowAdaptivity;
    bool                                    bVerbose;
    bool                                    bProfiler;
    bool                                    bVTK;
    Real                                    L;
    string                                  PatientFileName;
    int                                     rank;
    int                                     nprocs;
    Real                                    tumor_ic[3];

    
    static void _ic(Grid<W,B>& grid, int rank, string PatientFileName, Real& L, Real tumor_ic[3] );
    static void _icSphere3Parts(Grid<W,B>& grid, int rank, Real& L);
    double      _estimate_dt(double Diff_dt, double CFL);
    Real        _compute_maxvel();
    void        _computePressureSource(const int nParallelGranularity, const Real rho);
    void        _computeVelocities(BoundaryInfo* boundaryInfo ,const bool bMobility, std::vector<Real>* mobility);
    void        _reactionDiffusionStep(BoundaryInfo* boundaryInfo, const Real Dw, const Real Dg, const Real rho);
    void        _advectionConvectionStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, double dt);
    void		_dump(int counter);
    
public:
    Glioma_BrainDeformation(int argc, const char ** argv);
    ~Glioma_BrainDeformation();
    void run();
    
};


#endif /* defined(__Glioma_BrainDeformation__) */
