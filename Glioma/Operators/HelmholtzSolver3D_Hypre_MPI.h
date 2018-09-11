//
//  HelmholtzSolver3D_Hypre.h
//  GliomaSolverXcode
//
//  Created by Lipkova on 11/09/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
////  ------------------------------------
//   - ∇( ψ M∇u) + κψu = ψf
//     + periodic BC (autimatic in MRAG labs)
//  ------------------------------------
//     κ - relaxation factor
//     M  - mobility
//  ------------------------------------

// Assume MPI withine one node, or multiple nodes but whole system fits into memory of each node -> no message sendings supported yet


#ifndef GliomaSolverXcode_HelmholtzSolver4D_Hypre_h
#define GliomaSolverXcode_HelmholtzSolver4D_Hypre_h

#include "../Glioma_Types.h"

#include <_hypre_utilities.h>
#include <HYPRE_krylov.h>
#include <HYPRE_sstruct_ls.h>

class HelmholtzSolver3D_Hypre_MPI
{
    HYPRE_StructGrid     grid;
    HYPRE_StructStencil  stencil;
    HYPRE_StructMatrix   matrix;
    HYPRE_StructVector   rhs;
    HYPRE_StructVector   solution;
    HYPRE_StructSolver   solver;
    HYPRE_StructSolver   precond;
    
    int rank, size; // MPI_Commn_World
    
    Grid<W,B> * mrag_grid;
    
    bool bAlreadyAllocated;
    bool bVerbose;
    bool bCG;     // conjugated gradients with SGM are preconditioning
    bool bRelaxation;
    bool bMobility;
    
    Real kappaCSF, kappaWM, kappaGM;  // relaxation factors for diff. anatomies
    Real mWM, mGM, mCSF;      // hydralucit conductivuty (mobility)
    
    tbb::tick_count t1,t0;

    
#pragma mark Helpers
    void setup_hypre()
    {
        // deallocation
        if (bAlreadyAllocated)
            cleanUp();
        
        bAlreadyAllocated = true;
    }
    
    /* 1. Set up a grid. Each processor describes the piece of the grid called box it owns */
    void setup_structGrid()
    {
        // Set up grid
        HYPRE_StructGridCreate(MPI_COMM_WORLD, 3, &grid);
        
        printf("Hello from rank %d out of %d processes \n", rank, size);

        
    }
    
    
    void setup_structMatrix()
    {
        
    }
    
    void setup_Vectors()
    {
        
    }
    
    
    /* Set up a solver (See the Reference Manual for descriptions of all of the options.) */
    void setup_solver()
    {
        int n_pre   = 1;
        int n_post  = 1;
        
        if(bCG==0)
        {
            // use SMG solver
            HYPRE_StructSMGCreate(MPI_COMM_WORLD, &solver);
            HYPRE_StructSMGSetMemoryUse(solver, 0);
            HYPRE_StructSMGSetMaxIter(solver, 200);
            HYPRE_StructSMGSetTol(solver, 1.0e-12);
            HYPRE_StructSMGSetRelChange(solver, 0);
            HYPRE_StructSMGSetNonZeroGuess(solver);
            HYPRE_StructSMGSetNumPreRelax(solver, n_pre);
            HYPRE_StructSMGSetNumPostRelax(solver, n_post);
            
            /* Logging must be on to get iterations and residual norm info below */
            HYPRE_StructSMGSetPrintLevel(solver, 1);
            HYPRE_StructSMGSetLogging(solver, 1);
        }
        else
        {
            /* use CG with SMG as preconditioner */
            HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
            HYPRE_StructPCGSetMaxIter(solver, 200 );
            HYPRE_StructPCGSetTol(solver, 1.0e-12 );
            HYPRE_StructPCGSetTwoNorm(solver, 1 );
            HYPRE_StructPCGSetRelChange(solver, 0 );
            if(bVerbose) HYPRE_StructPCGSetPrintLevel(solver, 2 );
            
            HYPRE_StructSMGCreate(MPI_COMM_WORLD, &precond);
            HYPRE_StructSMGSetMemoryUse(precond, 0);
            HYPRE_StructSMGSetMaxIter(precond, 1);
            HYPRE_StructSMGSetTol(precond, 0.0);
            HYPRE_StructSMGSetZeroGuess(precond);
            HYPRE_StructSMGSetNumPreRelax(precond, n_pre);
            HYPRE_StructSMGSetNumPostRelax(precond, n_post);
            
            /* Logging must be on to get iterations and residual norm info below */
            HYPRE_StructSMGSetPrintLevel(precond, 0);
            HYPRE_StructSMGSetLogging(precond, 0);
            HYPRE_StructPCGSetPrecond(solver, HYPRE_StructSMGSolve,HYPRE_StructSMGSetup,precond);
        }
    }
    
    void solve()
    {
        
    }
    
    
    void getResultsOMP()
    {
        
    }
    
    void inline cleanUp()
    {
        HYPRE_StructGridDestroy(grid);
        HYPRE_StructStencilDestroy(stencil);
        HYPRE_StructMatrixDestroy(matrix);
        HYPRE_StructVectorDestroy(rhs);
        HYPRE_StructVectorDestroy(solution);
        
        
        if(bCG==0)
            HYPRE_StructSMGDestroy(solver);
        else
        {
            HYPRE_StructPCGDestroy(solver);
            HYPRE_StructSMGDestroy(precond);
        }
        
        bAlreadyAllocated = false;
        
    }
    
    
public:
    HelmholtzSolver3D_Hypre_MPI(int argc, const char ** argv): bAlreadyAllocated(false)
    {  }
    
    ~HelmholtzSolver3D_Hypre_MPI()
    { }
    
    void operator()(Grid<W,B>& input_grid, bool bVerbose=false, bool bCG=false, bool bRelaxation=false, std::vector<Real>* kappa = NULL, bool bMobility=false, std::vector<Real>* mobility = NULL)
    {
        mrag_grid           = &input_grid;
        this->bVerbose      = bVerbose;
        this->bCG           = bCG;
        this->bRelaxation   = bRelaxation;
        
        
        if (bRelaxation)
            assert (kappa!=NULL);
        
        this->kappaCSF      =  (bRelaxation) ? (*kappa)[0] : 1. ;
        this->kappaWM       =  (bRelaxation) ? (*kappa)[1] : 1. ;
        this->kappaGM       =  (bRelaxation) ? (*kappa)[2] : 1. ;
        
        
        if(bMobility)
            assert (mobility!=NULL);
        
        this-> mCSF         = (bMobility) ? (*mobility)[0] : 1. ;
        this-> mWM          = (bMobility) ? (*mobility)[1] : 1. ;
        this-> mGM          = (bMobility) ? (*mobility)[2] : 1. ;
        
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        
        setup_hypre();
        setup_structGrid();
        setup_structMatrix();
        setup_Vectors();
        setup_solver();
        solve();
        getResultsOMP();
        cleanUp();
        
    }
    
};


#endif
