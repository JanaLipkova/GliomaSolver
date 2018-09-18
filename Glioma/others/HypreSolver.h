//
//  HypreSolver.h
//  GliamoSolverXcode_kraken
//
//  Created by Lipkova on 18/09/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

#pragma once
#include "../Glioma_Types.h"
#include "../Glioma.h"
#include "../../MRAG/MRAGcore/MRAGrid.h"

class HypreSolver
{
public:
    virtual void setup_hypre( Grid<W,B>& input_grid, int rank, int nprocs, bool bVerbose=false, bool bCG=false, bool bRelaxation=false, std::vector<Real>* kappa = NULL, bool bMobility=false, std::vector<Real>* mobility = NULL  ) = 0;
    virtual void solve() = 0;
    virtual void clean() = 0;
    
    virtual ~HypreSolver(){}
};
