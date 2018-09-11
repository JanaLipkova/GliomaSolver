//
//  HelmholtzSolver3D_Hypre_MPI.h
//  GliomaSolver
//
//  Created by Lipkova on 11/09/18.
//  based on Hypre example codes
//  Copyright (c) 2018 Lipkova. All rights reserved.
//  ------------------------------------
//   - ∇( ψ ∇u) + ψu = ψf
//     + periodic BC (automatic in MRAG)
//  ------------------------------------
//
// Assume 1) uniform grid (no refinemnet)
//        2) mod(#MRAG_blocks,MPI_process)==0


#include <mpi.h>
#include <_hypre_utilities.h>
#include <HYPRE_krylov.h>
#include <HYPRE_sstruct_ls.h>
#include "../Glioma_Types.h"


class HelmholtzSolver3D_Hypre_MPI
{
    HYPRE_StructGrid     grid;
    HYPRE_StructStencil  stencil;
    HYPRE_StructMatrix   matrix;
    HYPRE_StructVector   rhs;
    HYPRE_StructVector   solution;
    HYPRE_StructSolver   solver;
    HYPRE_StructSolver   precond;
    
    bool bAlreadyAllocated;
    bool bVerbose;
    bool bRelaxation;
    bool bMobility;
    
    Real kappaCSF, kappaWM, kappaGM;  // relaxation factors for diff. anatomies
    Real mWM, mGM, mCSF;      // hydralucit conductivuty (mobility)
    
    tbb::tick_count t1,t0;
    
    Grid<W,B> * mrag_grid;
    
    int myid, num_procs;
 
    /* Set up a solver (See the Reference Manual for descriptions of all of the options.) */
    void _setup_solver()
    {
        int n_pre   = 1;
        int n_post  = 1;
        
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
    
    void inline _cleanUp()
    {
        HYPRE_StructGridDestroy(grid);
        HYPRE_StructStencilDestroy(stencil);
        HYPRE_StructMatrixDestroy(matrix);
        HYPRE_StructVectorDestroy(rhs);
        HYPRE_StructVectorDestroy(solution);
        HYPRE_StructPCGDestroy(solver);
        HYPRE_StructSMGDestroy(precond);
        
        bAlreadyAllocated = false;
    }
    
public:
    HelmholtzSolver3D_Hypre_MPI(int argc, const char ** argv): bAlreadyAllocated(false)
    { }
    
    ~HelmholtzSolver3D_Hypre_MPI()
    { }
    
    void _setup_hypre()
    {
        if (bAlreadyAllocated)
            _cleanUp();
        
        bAlreadyAllocated = true;
        if (bVerbose) printf("done with setup!\n"); //exit(0);
    }
    
    void operator()(Grid<W,B>& input_grid, bool bVerbose=false, bool bRelaxation=false, std::vector<Real>* kappa = NULL, bool bMobility=false, std::vector<Real>* mobility = NULL)
    {
        mrag_grid           = &input_grid;
        this->bVerbose      = bVerbose;
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


        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
        printf("Hello world from processor rank %d out of %d processors\n", myid, num_procs);
        
        
        printf("Calling set-up Hypre \n");
        _setup_hypre();
        
        printf("Calling set-up Solver \n");
        _setup_solver();
        
        printf("Calling clean-up \n");
        _cleanUp();
      
    }
    
};



