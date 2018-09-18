//
//  HelmholtzSolver2D_Hypre.h
//  GliomaXcode
//
//  Created by Lipkova on 15/04/15,
//  based on Christian Conti code and Hypre example codes
//  Copyright (c) 2015 Lipkova. All rights reserved.
//
//  ------------------------------------
//   - ∇( ψ M∇u) + κψu = ψf
//     + periodic BC
//     κ - relaxation factor
//     M  - mobility
//  ------------------------------------


#ifndef HelmholtzSolver2D_Hypre_h
#define HelmholtzSolver2D_Hypre_h

#include <_hypre_utilities.h>
#include <HYPRE_krylov.h>
#include <HYPRE_sstruct_ls.h>

class HelmholtzSolver2D_Hypre
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
    bool bCG;              // use for conjugated gradients with SGM are preconditioning
    bool bRelaxation;
    bool bMobility;
    
    Real kappaWM, kappaGM, kappaCSF;  // relaxation factors for diff. anatomies
    Real mWM, mGM, mCSF;      // hydralucit conductivuty (mobility)
    
    int GridsizeX, GridsizeY;
    Grid<W,B> * mrag_grid;
    
    /* Set up a Struct Matrix */
    inline void _setupMatrix()
    {
        /* Create an empty matrix object */
        HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &matrix);
        
        /* Indicate that the matrix coefficients are ready to be set */
        HYPRE_StructMatrixInitialize(matrix);

        /* Set the matrix coefficients.  Each processor assigns coefficients for the boxes in the grid that it owns. Note that the coefficients associated with each stencil entry may vary from grid point to grid point if desired. */
        int ilower[2]           = { 0, 0 };
        int iupper[2]           = { GridsizeX-1, GridsizeY-1 };
        int stencil_indices[5]  = {0,1,2,3,4}; // labels for the stencil entries compatible with offsets defined in _hypre_setup()
        
        const int nentries      = 5;
        const int size2         = GridsizeX * GridsizeY;
        int nvalues             = size2 * nentries;        // We have size2 grid points, each with 5 stencil entries
        double * values         = new double [nvalues];
        
        _fillMatrixOMP<BlockLab<B> >(values);
        
        HYPRE_StructMatrixSetBoxValues(matrix, ilower, iupper, nentries, stencil_indices, values);
        
        /* This is a collective call finalizing the matrix assembly.The matrix is now ``ready to be used'' */
        HYPRE_StructMatrixAssemble(matrix);
        
        delete [] values;
    }
    
    template<typename BlockLabType>
    void _fillMatrixOMP(double * values)
    {
        vector<BlockInfo> vInfo = mrag_grid->getBlocksInfo();
        const BlockCollection<B>& coll = mrag_grid->getBlockCollection();
        const BoundaryInfo& binfo=mrag_grid->getBoundaryInfo();
        
        const int stencilStart[3] = {-1, -1, 0};
        const int stencilEnd[3]   = {+2, +2, 1};
        
        
#pragma omp parallel
        {
            BlockLabType lab;
            lab.prepare(coll, binfo, stencilStart, stencilEnd);
            
            
#pragma omp for schedule(static)
            for(int i=0;i<vInfo.size();i++)
            {
                lab.load(vInfo[i]);
                const BlockInfo info = vInfo[i];
                B& block = mrag_grid->getBlockCollection()[info.blockID];
                
                double h       = info.h[0];
                double h2      = h*h;
                
                Real m, mS, mN, mE, mW;
                
                for(int iy=0; iy<B::sizeY; ++iy)
                    for(int ix=0; ix<B::sizeX; ++ix)
                    {
                        // get global index
                        const int gix = ix + info.index[0] * B::sizeX;
                        const int giy = iy + info.index[1] * B::sizeY;
                        const int idx = gix*5 + giy * 5 * blocksPerDimension*B::sizeX;
                        
                        
                        // compute mobility components
                        if( (mWM == mGM) && (mWM == mCSF) )
                            m = mS = mN = mW = mE = mCSF;
                        else
                        {
                            mS  = mWM * lab(ix,  iy-1).p_w  + mGM * lab(ix,  iy-1).p_g  + mCSF * lab(ix,  iy-1).p_csf ;
                            mN  = mWM * lab(ix,  iy+1).p_w  + mGM * lab(ix,  iy+1).p_g  + mCSF * lab(ix,  iy+1).p_csf ;
                            mW  = mWM * lab(ix-1,iy  ).p_w  + mGM * lab(ix-1,iy  ).p_g  + mCSF * lab(ix-1,iy  ).p_csf ;
                            mE  = mWM * lab(ix+1,iy  ).p_w  + mGM * lab(ix+1,iy  ).p_g  + mCSF * lab(ix+1,iy  ).p_csf ;
                            m   = mWM * lab(ix,  iy  ).p_w  + mGM * lab(ix,  iy  ).p_g  + mCSF * lab(ix,  iy  ).p_csf ;
                        }

                        // get entries of matrix
                        Real pffS = lab(ix,  iy-1).pff * mS;
                        Real pffN = lab(ix,  iy+1).pff * mN;
                        Real pffW = lab(ix-1,iy  ).pff * mW;
                        Real pffE = lab(ix+1,iy  ).pff * mE;
                        Real pff  = lab(ix,iy).pff * m;

                        
                        // approximate intermidiet points
                        _mean(pff, pffW, pffE, pffS, pffN);
                        
                        Real kappa = kappaWM * lab(ix,iy).p_w  + kappaGM * lab(ix,iy).p_g  + kappaCSF * lab(ix,iy).p_csf;
                        
                        // fill in vector of matrix values
                        values[idx  ] =   pffW + pffE + pffS + pffN + kappa *pff*h2;
                        values[idx+1] = - pffW;
                        values[idx+2] = - pffE;
                        values[idx+3] = - pffS;
                        values[idx+4] = - pffN;
                        
                    }
            }
        }
        
    }
    
    inline void _mean(Real& psi, Real& psiW, Real& psiE, Real& psiS, Real& psiN)
    {
        psiW = 0.5 * (psi + psiW);
        psiE = 0.5 * (psi + psiE);
        psiS = 0.5 * (psi + psiS);
        psiN = 0.5 * (psi + psiN);
    }
    
    inline void _harmonicAvg(Real& psi, Real& psiW, Real& psiE, Real& psiS, Real& psiN)
    {
        psiW = 2. * psi * psiW / (psi + psiW);
        psiE = 2. * psi * psiE / (psi + psiE);
        psiS = 2. * psi * psiS / (psi + psiS);
        psiN = 2. * psi * psiN / (psi + psiN);
    }

    /* Set up Struct Vectors for rhs and solution */
    inline void _setupVectors()
    {
        /* Create an empty vector object */
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &rhs);
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &solution);
        
        /* Indicate that the vector coefficients are ready to be set */
        HYPRE_StructVectorInitialize(rhs);
        HYPRE_StructVectorInitialize(solution);
        
        /* Set the vector coefficients */
        int ilower[2]       = { 0, 0 };
        int iupper[2]       = { GridsizeX-1, GridsizeY-1 };
        const int size2     = GridsizeX * GridsizeY;
        
        double * valuesRhs  = new double[size2];
        double * valuesSol  = new double[size2];

        
        /* fill right hand side */
        vector<BlockInfo> vInfo = mrag_grid->getBlocksInfo();
        
#pragma omp parallel for schedule(static)
        for(int i=0; i<vInfo.size(); i++)
        {
            const BlockInfo info = vInfo[i];
            B& block = mrag_grid->getBlockCollection()[info.blockID];
            
            double h       = info.h[0];
            double h2      = h*h;
            
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    const int gix = ix + info.index[0] * B::sizeX;
                    const int giy = iy + info.index[1] * B::sizeY;
                    const int idx = gix + giy * blocksPerDimension*B::sizeX;
                    
                    valuesRhs[idx] = h2 * block(ix,iy).f;
                    valuesSol[idx] =      block(ix,iy).p;

                }
        }
        
        HYPRE_StructVectorSetBoxValues(rhs, ilower, iupper, valuesRhs);
        HYPRE_StructVectorSetBoxValues(solution, ilower, iupper, valuesSol);
        
        /* This is a collective call finalizing the vector assembly.
         The vectors are now ``ready to be used'' */
        HYPRE_StructVectorAssemble(rhs);
        HYPRE_StructVectorAssemble(solution);
        
        delete [] valuesRhs;
        delete [] valuesSol;
    }
    
    /* Set up a solver (See the Reference Manual for descriptions of all of the options.) */
    void _setup_solver()
    {
        int n_pre   = 1;
        int n_post  = 1;
        
        if(bCG==0)
        {
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
            /* use symmetric SMG as preconditioner */
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
            HYPRE_StructSMGSetNonZeroGuess(precond);
            HYPRE_StructSMGSetNumPreRelax(precond, n_pre);
            HYPRE_StructSMGSetNumPostRelax(precond, n_post);
            HYPRE_StructSMGSetPrintLevel(precond, 0);
            HYPRE_StructSMGSetLogging(precond, 0);
            HYPRE_StructPCGSetPrecond(solver, HYPRE_StructSMGSolve,HYPRE_StructSMGSetup,precond);
        }
    }
    
    void _solveSystem()
    {
        int num_iterations;
        double final_res_norm;
        
        /* Setup and solve */
        if(bCG==0)
        {
            HYPRE_StructSMGSetup(solver, matrix, rhs, solution);
            HYPRE_StructSMGSolve(solver, matrix, rhs, solution);
            
            /* Get info */
            HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
            HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
        }
        else
        {
            HYPRE_StructPCGSetup(solver, matrix, rhs, solution );
            HYPRE_StructPCGSolve(solver, matrix, rhs, solution);
            
            HYPRE_StructPCGGetNumIterations( solver, &num_iterations );
            HYPRE_StructPCGGetFinalRelativeResidualNorm( solver, &final_res_norm );
        }
        
        if(bVerbose)
        {
            printf("=========== HYPRE ===========\n");
            printf("Iterations = %d\n", num_iterations);
            printf("Final Relative Residual Norm = %g\n", final_res_norm);
            printf("============================\n");
        }
    }
    
    void _getResultsOMP()
    {
        // get solution
        int ilower[2]           = { 0, 0 };
        int iupper[2]           = { GridsizeX-1, GridsizeY-1 };
        const int size2         = GridsizeX * GridsizeY;
        double values[size2];
        HYPRE_StructVectorGetBoxValues(solution, ilower, iupper, values);
        
        // write solution to the grid
        vector<BlockInfo> vInfo = mrag_grid->getBlocksInfo();
        
        
#pragma omp parallel for schedule(static)
        for(int i=0; i<vInfo.size(); i++)
        {
            const BlockInfo info = vInfo[i];
            B& block = mrag_grid->getBlockCollection()[info.blockID];
                        
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    const int gix = ix + info.index[0] * B::sizeX;
                    const int giy = iy + info.index[1] * B::sizeY;
                    const int idx = gix + giy * blocksPerDimension*B::sizeX;
                    
                    block(ix,iy).p = values[idx];
                }
        }
        
    }
    
    
public:
    HelmholtzSolver2D_Hypre(): bAlreadyAllocated(false)
    {  }
    
    ~HelmholtzSolver2D_Hypre()
    { }
    
    void setup_hypre(Grid<W,B>& input_grid, int rank, int nprocs, bool bVerbose=false, bool bCG=false, bool bRelaxation=false, std::vector<Real>* kappa = NULL, bool bMobility=false, std::vector<Real>* mobility = NULL)
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
        
        //0. deallocation
        if (bAlreadyAllocated)
            clean();
        
        GridsizeX = blocksPerDimension * B::sizeX;
        GridsizeY = blocksPerDimension * B::sizeY;
        assert(B::sizeZ == 1);
        
        /* 1. Set up a grid. Each processor describes the piece of the grid it owns */
        HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);
        
        int ilower[2] = { 0, 0 };
        int iupper[2] = { GridsizeX-1, GridsizeY-1 };
        HYPRE_StructGridSetExtents(grid, ilower, iupper);
        
        int periodicity[2] = { GridsizeX,GridsizeY };
        HYPRE_StructGridSetPeriodic(grid, periodicity);
        
        HYPRE_StructGridAssemble(grid);
        
        /* 2. Define the discretization stencil */
        HYPRE_StructStencilCreate(2, 5, &stencil);
        
        int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};
        
        for (int entry = 0; entry < 5; entry++)
            HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
        
        /*3. Set up solver */
        _setup_solver();
        
        bAlreadyAllocated = true;
        if (bVerbose) printf("done with setup!\n"); //exit(0);
    }
    
    void inline solve()
    {
        _setupMatrix();
        _setupVectors();
        _solveSystem();
        _getResultsOMP();
    }
    
    void inline clean()
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
    
};


#endif
