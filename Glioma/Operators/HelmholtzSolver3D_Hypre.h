//
//  HelmholtzSolver3D_Hypre.h
//  GliomaXcode
//
//  Created by Lipkova on 23/04/15.
//  based on Christian Conti code and Hypre example codes
//  Copyright (c) 2015 Lipkova. All rights reserved.
//  ------------------------------------
//   - ∇( ψ ∇u) + ψu = ψf
//     + periodic BC (automatic in MRAG)
//  ------------------------------------

#ifndef HelmholtzSolver3D_Hypre_h
#define HelmholtzSolver3D_Hypre_h

#include "../Glioma_Types.h"

#include <_hypre_utilities.h>
#include <HYPRE_krylov.h>
#include <HYPRE_sstruct_ls.h>

class HelmholtzSolver3D_Hypre
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
    bool bCG;     // conjugated gradients with SGM are preconditioning
    bool bRelaxation;
    bool bMobility;

    
    Real kappaCSF, kappaWM, kappaGM;  // relaxation factors for diff. anatomies
    Real mWM, mGM, mCSF;      // hydralucit conductivuty (mobility)

    tbb::tick_count t1,t0;
    
    int GridsizeX, GridsizeY,GridsizeZ ;
    Grid<W,B> * mrag_grid;
    
    /* Set up a Struct Matrix */
    inline void _setupMatrix()
    {
        /* Create an empty matrix object */
        HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &matrix);
        
        /* Indicate that the matrix coefficients are ready to be set */
        HYPRE_StructMatrixInitialize(matrix);
        
        /* Set the matrix coefficients.  Each processor assigns coefficients for the boxes in the grid that it owns. Note that the coefficients associated with each stencil entry may vary from grid point to grid point if desired. */
        int ilower[3]           = { 0, 0, 0 };
        int iupper[3]           = {GridsizeX-1, GridsizeY-1, GridsizeZ-1 };
        int stencil_indices[7]  = {0,1,2,3,4,5,6}; // labels for the stencil entries, compatible with offsets defined in _hypre_setup()
        
        const int nentries      = 7;
        const int size3         = GridsizeX * GridsizeY * GridsizeZ;
        const int nvalues       = size3 * nentries;
        double * values         = new double [nvalues];
        
        _fillMatrixOMP<BlockLab<B> >(values);
        
        HYPRE_StructMatrixSetBoxValues(matrix, ilower, iupper, nentries, stencil_indices, values);
        
        /* This is a collective call finalizing the matrix assembly.
         The matrix is now ``ready to be used'' */
        HYPRE_StructMatrixAssemble(matrix);
        
        delete [] values;
    }
    
    template<typename BlockLabType>
    void _fillMatrixOMP(double * values)
    {
        vector<BlockInfo> vInfo = mrag_grid->getBlocksInfo();
        const BlockCollection<B>& coll = mrag_grid->getBlockCollection();
        const BoundaryInfo& binfo=mrag_grid->getBoundaryInfo();
        
        const int stencilStart[3] = { -1, -1, -1};
        const int stencilEnd[3]   = { +2, +2, +2};
        
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
                Real m, mS, mN, mW, mE, mF, mB;
                
                for(int iz=0; iz<B::sizeZ; ++iz)
                    for(int iy=0; iy<B::sizeY; ++iy)
                        for(int ix=0; ix<B::sizeX; ++ix)
                        {
                            // get global index
                            const int gix = ix + info.index[0] * B::sizeX;
                            const int giy = iy + info.index[1] * B::sizeY;
                            const int giz = iz + info.index[2] * B::sizeZ;
                            
                            const int idx = gix*7 + giy * 7 * blocksPerDimension * B::sizeX
                            + giz * 7 * blocksPerDimension * B::sizeX * blocksPerDimension*B::sizeY;
                            
                            
                            // compute mobility components
                            if ((mWM == mGM) && (mWM == mCSF))
                                m = mS = mN = mW = mE = mF = mB = mCSF;
                            else
                            {
                                mB  = mWM * lab(ix,  iy,  iz-1).p_w + mGM * lab(ix,  iy,  iz-1).p_g + mCSF * lab(ix,  iy,  iz-1).p_csf;
                                mF  = mWM * lab(ix,  iy,  iz+1).p_w + mGM * lab(ix,  iy,  iz+1).p_g + mCSF * lab(ix,  iy,  iz+1).p_csf;
                                mS  = mWM * lab(ix,  iy-1,iz  ).p_w + mGM * lab(ix,  iy-1,iz  ).p_g + mCSF * lab(ix,  iy-1,iz  ).p_csf;
                                mN  = mWM * lab(ix,  iy+1,iz  ).p_w + mGM * lab(ix,  iy+1,iz  ).p_g + mCSF * lab(ix,  iy+1,iz  ).p_csf;
                                mW  = mWM * lab(ix-1,iy  ,iz  ).p_w + mGM * lab(ix-1,iy  ,iz  ).p_g + mCSF * lab(ix-1,iy  ,iz  ).p_csf;
                                mE  = mWM * lab(ix+1,iy  ,iz  ).p_w + mGM * lab(ix+1,iy  ,iz  ).p_g + mCSF * lab(ix+1,iy  ,iz  ).p_csf;
                                m   = mWM * lab(ix,  iy,  iz  ).p_w + mGM * lab(ix,  iy,  iz  ).p_g + mCSF * lab(ix,  iy,  iz  ).p_csf;
                            }
                            
                            // get entries of matrix
                            Real pff  = lab(ix,iy,iz).pff * m;
                            Real pffB = lab(ix  ,iy  ,iz-1).pff * mB;  //back
                            Real pffF = lab(ix  ,iy  ,iz+1).pff * mF;  //front
                            Real pffS = lab(ix,  iy-1,iz  ).pff * mS;
                            Real pffN = lab(ix,  iy+1,iz  ).pff * mN;
                            Real pffW = lab(ix-1,iy  ,iz  ).pff * mW;
                            Real pffE = lab(ix+1,iy  ,iz  ).pff * mE;
                            
                            // approximate intermidiet points
                            _mean(pff, pffW, pffE, pffS, pffN, pffB, pffF);
                            
                            Real kappa = kappaWM * lab(ix,iy,iz).p_w + kappaGM * lab(ix,iy,iz).p_g  + kappaCSF * lab(ix,iy,iz).p_csf;
                            
                            // fill in vector of matrix values
                            values[idx  ] =   pffW + pffE + pffS + pffN + pffB + pffF + kappa*pff*h2;
                            values[idx+1] = - pffW;
                            values[idx+2] = - pffE;
                            values[idx+3] = - pffS;
                            values[idx+4] = - pffN;
                            values[idx+5] = - pffB;
                            values[idx+6] = - pffF;
                            
                        }
            }
        }
        
    }
    
    inline void _mean(Real& psi, Real& psiW, Real& psiE, Real& psiS, Real& psiN, Real& psiB, Real& psiF)
    {
        psiW = 0.5 * (psi + psiW);
        psiE = 0.5 * (psi + psiE);
        psiS = 0.5 * (psi + psiS);
        psiN = 0.5 * (psi + psiN);
        psiB = 0.5 * (psi + psiB);
        psiF = 0.5 * (psi + psiF);
        
    }
    
    inline void _harmonicAvg(Real& psi, Real& psiW, Real& psiE, Real& psiS, Real& psiN, Real& psiB, Real& psiF)
    {
        psiW = 2. * psi * psiW / (psi + psiW);
        psiE = 2. * psi * psiE / (psi + psiE);
        psiS = 2. * psi * psiS / (psi + psiS);
        psiN = 2. * psi * psiN / (psi + psiN);
        psiB = 2. * psi * psiB / (psi + psiB);
        psiF = 2. * psi * psiF / (psi + psiF);
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
        int ilower[3]       = { 0, 0,0 };
        int iupper[3]       = { GridsizeX-1, GridsizeY-1, GridsizeZ-1 };
        const int size3     = GridsizeX * GridsizeY * GridsizeZ;
        
        double * valuesRhs     = new double[size3];
        double * valuesSol     = new double[size3];
        
        vector<BlockInfo> vInfo = mrag_grid->getBlocksInfo();
        
#pragma omp parallel for schedule(static)
        for(int i=0; i<vInfo.size(); i++)
        {
            const BlockInfo info = vInfo[i];
            B& block = mrag_grid->getBlockCollection()[info.blockID];
            
            double h       = info.h[0];
            double h2      = h*h;
            
            for(int iz=0; iz<B::sizeZ; iz++)
                for(int iy=0; iy<B::sizeY; iy++)
                    for(int ix=0; ix<B::sizeX; ix++)
                    {
                        const int gix = ix + info.index[0] * B::sizeX;
                        const int giy = iy + info.index[1] * B::sizeY;
                        const int giz = iz + info.index[2] * B::sizeZ;
                        
                        const int idx = gix + giy * blocksPerDimension*B::sizeX
                        + giz * blocksPerDimension*B::sizeX * blocksPerDimension*B::sizeY;
                        
                        valuesRhs[idx] = h2 * block(ix,iy,iz).f;
                        valuesSol[idx] =      block(ix,iy,iz).p;
                        
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
    
    void _solve()
    {
        int num_iterations;
        double final_res_norm;
        
        /* Setup and solve */
        if(bCG==0)
        {
            HYPRE_StructSMGSetup(solver, matrix, rhs, solution);
            HYPRE_StructSMGSolve(solver, matrix, rhs, solution);
            
            HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
            HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
        }
        else
        {
            
            HYPRE_StructPCGSetup(solver, matrix, rhs, solution );
            HYPRE_StructPCGSolve(solver, matrix, rhs, solution);
            
            /* Get info */
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
        int ilower[3]       = { 0, 0, 0 };
        int iupper[3]       = { GridsizeX-1, GridsizeY-1, GridsizeZ-1 };
        const int size3     = GridsizeX * GridsizeY * GridsizeZ;
        double * values     = new double[size3];
        HYPRE_StructVectorGetBoxValues(solution, ilower, iupper, values);
        
        vector<BlockInfo> vInfo = mrag_grid->getBlocksInfo();
        
#pragma omp parallel for schedule(static)
        for(int i=0; i<vInfo.size(); i++)
        {
            const BlockInfo info = vInfo[i];
            B& block = mrag_grid->getBlockCollection()[info.blockID];
            
            for(int iz=0; iz<B::sizeZ; iz++)
                for(int iy=0; iy<B::sizeY; iy++)
                    for(int ix=0; ix<B::sizeX; ix++)
                    {
                        const int gix = ix + info.index[0] * B::sizeX;
                        const int giy = iy + info.index[1] * B::sizeY;
                        const int giz = iz + info.index[2] * B::sizeZ;
                        
                        const int idx = gix + giy * blocksPerDimension * B::sizeX
                        + giz * blocksPerDimension * B::sizeX * blocksPerDimension*B::sizeY;
                        
                        block(ix,iy,iz).p = values[idx];
                    }
        }
        
        delete [ ] values;
    }
    
    void inline _cleanUp()
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
    HelmholtzSolver3D_Hypre(int argc, const char ** argv): bAlreadyAllocated(false)
    {  }
    
    ~HelmholtzSolver3D_Hypre()
    { }
    
    void setup_hypre()
    {
        //0. deallocation
        //1. setup della grid di merda
        //2. setup dello stencil del cazzo
        
        //0.
        if (bAlreadyAllocated)
            _cleanUp();
        
        GridsizeX = blocksPerDimension * B::sizeX;
        GridsizeY = blocksPerDimension * B::sizeY;
        GridsizeZ = blocksPerDimension * B::sizeZ;
        
        /* 1. Set up a grid. Each processor describes the piece of the grid it owns */
        HYPRE_StructGridCreate(MPI_COMM_WORLD, 3, &grid);
        
        int ilower[3] = { 0, 0, 0 };
        int iupper[3] = { GridsizeX-1, GridsizeY-1, GridsizeZ -1 };
        HYPRE_StructGridSetExtents(grid, ilower, iupper);
        
        int periodicity[3] = { GridsizeX,GridsizeY, GridsizeZ };
        HYPRE_StructGridSetPeriodic(grid, periodicity);
        
        HYPRE_StructGridAssemble(grid);
        
        /* 2. Define the discretization stencil */
        HYPRE_StructStencilCreate(3, 7, &stencil);
        
        int offsets[7][3] = {{0,0,0}, {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1},{0,0,1}};
        
        for (int entry = 0; entry < 7; entry++)
            HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
        
        bAlreadyAllocated = true;
        if (bVerbose) printf("done with setup!\n"); //exit(0);
    }
    
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

        
        if(bVerbose) printf("Relaxation factors: kappaCSF=%f  kappaWM=%f, kappaGM=%f\n", kappaCSF, kappaWM, kappaGM);
        if(bVerbose) printf("Mobility factors:  mWM=%f, mGM=%f, mCSF=%f\n", mWM, mGM, mCSF);

        
        setup_hypre();
        _setupMatrix();
        _setupVectors();
        _setup_solver();
        _solve();
        _getResultsOMP();
        
        _cleanUp();

        
    }
    
};

#endif
