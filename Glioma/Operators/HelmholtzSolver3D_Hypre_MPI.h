//
//  HelmholtzSolver3D_Hypre_MPI.h
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


#ifndef GliomaSolverXcode_HelmholtzSolver3D_Hypre_MPI_h
#define GliomaSolverXcode_HelmholtzSolver3D_Hypre_MPI_h
class HelmholtzSolver3D_Hypre_MPI
{
    HYPRE_StructGrid     hypre_grid;
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
    
    Grid<W,B> * mrag_grid;
    
    int GridsizeX, GridsizeY,GridsizeZ ;
    int nMRAGblocks;
    
    
    int blocksPerProcesor;
    int rank, nprocs;
    
//    inline void _setup_grid()
//    {
//        /* 1. Set up a grid. Each processor describes the piece of the grid it owns */
//        HYPRE_StructGridCreate(MPI_COMM_WORLD, 3, &grid);
//        
//        int ilower[3] = { 0, 0, 0 };
//        int iupper[3] = { GridsizeX-1, GridsizeY-1, GridsizeZ -1 };
//        HYPRE_StructGridSetExtents(grid, ilower, iupper);
//        
//        int periodicity[3] = { GridsizeX,GridsizeY, GridsizeZ };
//        HYPRE_StructGridSetPeriodic(grid, periodicity);
//        
//        HYPRE_StructGridAssemble(grid);
//
//    }
    
    /* Set up the grid */
    inline void _setup_grid()
    {
        // 1) Creates grid object
        HYPRE_StructGridCreate(MPI_COMM_WORLD, 3, &hypre_grid);
        
        // 2) For each processor set the mrag blocks it will own
        vector<BlockInfo> vInfo = mrag_grid->getBlocksInfo();
        
        for(int i = 0; i < blocksPerProcesor; i++)
        {
            int blockID     = rank + (blocksPerProcesor - 1) * rank + i;
            BlockInfo& info = vInfo[blockID];
            B& block        = mrag_grid->getBlockCollection()[info.blockID];
            
            int ilower[3]   = {info.index[0] * B::sizeX,
                info.index[1] * B::sizeY,
                info.index[2] * B::sizeZ };
            
            int iupper[3]   = { ilower[0] + B::sizeX - 1,
                ilower[1] + B::sizeY - 1,
                ilower[2] + B::sizeZ - 1};
            
            HYPRE_StructGridSetExtents(hypre_grid, ilower, iupper);
        }
        
        // 3) Set periodicity:
        int periodicity[3] = {GridsizeX,GridsizeY, GridsizeZ };
        HYPRE_StructGridSetPeriodic(hypre_grid, periodicity);
        
        // 4) Assemble grid (collective call)
        HYPRE_StructGridAssemble(hypre_grid);
    }
    
    inline void _setup_stencil()
    {
        HYPRE_StructStencilCreate(3, 7, &stencil);
        int offsets[7][3] = {{0,0,0}, {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1},{0,0,1}};
        
        for (int entry = 0; entry < 7; entry++)
            HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
    }
    
    
    /* Set up a Struct Matrix */
    inline void _setupMatrix()
    {
        /* Create an empty matrix object */
        HYPRE_StructMatrixCreate(MPI_COMM_WORLD, hypre_grid, stencil, &matrix);
        
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
    
    
//    /* Set up a Struct Matrix */
//    inline void _setupMatrix()
//    {
//        /* Create an empty matrix object + indicate it is ready to be set*/
//        HYPRE_StructMatrixCreate(MPI_COMM_WORLD, hypre_grid, stencil, &matrix);
//        HYPRE_StructMatrixInitialize(matrix);
//        
//        /* For each processor set the matrix part based on the MRAG blocks it owns*/
//        
//        int stencil_indices[7]  = {0,1,2,3,4,5,6}; // labels for the stencil entries, compatible with offsets defined in _hypre_setup()
//        
//        const int nentries      = 7;
//        const int blockSize     = B::sizeX * B::sizeY * B::sizeZ;
//        const int nvalues       = nentries * blockSize;
//
//        vector<BlockInfo> vInfo = mrag_grid->getBlocksInfo();
//        const BlockCollection<B>& coll = mrag_grid->getBlockCollection();
//        const BoundaryInfo& binfo=mrag_grid->getBoundaryInfo();
//        
//        const int stencilStart[3] = { -1, -1, -1};
//        const int stencilEnd[3]   = { +2, +2, +2};
//        
//
//        for(int i = 0; i < blocksPerProcesor; i++)
//        {
//            int bID     = rank + (blocksPerProcesor - 1) * rank + i;
//            BlockInfo& info = vInfo[bID];
//            B& block        = mrag_grid->getBlockCollection()[info.blockID];
//            
//            int ilower[3]   = {info.index[0] * B::sizeX,
//                info.index[1] * B::sizeY,
//                info.index[2] * B::sizeZ };
//            
//            int iupper[3]   = { ilower[0] + B::sizeX - 1,
//                ilower[1] + B::sizeY - 1,
//                ilower[2] + B::sizeZ - 1};
//            
//            vector<double> values(nvalues);
//
//            
//            BlockLab<B> lab;
//            lab.prepare(coll, binfo, stencilStart, stencilEnd);
//            
//            double h       = info.h[0];
//            double h2      = h*h;
//            Real m, mS, mN, mW, mE, mF, mB;
//            
//            for(int iz=0; iz<B::sizeZ; ++iz)
//                for(int iy=0; iy<B::sizeY; ++iy)
//                    for(int ix=0; ix<B::sizeX; ++ix)
//                    {
//                        // get global index
//                        const int gix = ix + info.index[0] * B::sizeX;
//                        const int giy = iy + info.index[1] * B::sizeY;
//                        const int giz = iz + info.index[2] * B::sizeZ;
//                        
//                        const int idx = gix*7 + giy * 7 * blocksPerDimension * B::sizeX
//                        + giz * 7 * blocksPerDimension * B::sizeX * blocksPerDimension*B::sizeY;
//                        
//                        values[idx  ] =   4;
//                        values[idx+1] = - 1;
//                        values[idx+2] = - 1;
//                        values[idx+3] = - 1;
//                        values[idx+4] = - 1;
//                        values[idx+5] = - 1;
//                        values[idx+6] = - 1;
//                    }
//
//            HYPRE_StructMatrixSetBoxValues(matrix, ilower, iupper, nentries, stencil_indices, &values[0]);
//
//        }
//        
//        /* This is a collective call finalizing the matrix assembly.
//         The matrix is now ``ready to be used'' */
//        HYPRE_StructMatrixAssemble(matrix);
//    }
//    
////    template<typename BlockLabType>
//    //    void _fillMatrix(double * values, int blockID)
//    template<typename BlockLabType>
//    void _fillMatrix(std::vector<double> & values, int bID )
//    {
//        vector<BlockInfo> vInfo = mrag_grid->getBlocksInfo();
//        const BlockCollection<B>& coll = mrag_grid->getBlockCollection();
//        const BoundaryInfo& binfo=mrag_grid->getBoundaryInfo();
//        
//        const int stencilStart[3] = { -1, -1, -1};
//        const int stencilEnd[3]   = { +2, +2, +2};
//        
//        printf("Inside fillMatrix with %d blockID \n", bID);
//        
//        BlockLabType lab;
//        //            lab.prepare(coll, binfo, stencilStart, stencilEnd, false);
//        lab.prepare(coll, binfo, stencilStart, stencilEnd);
//        
//        lab.load(vInfo[bID]);
//        const BlockInfo info = vInfo[bID];
//        B& block = mrag_grid->getBlockCollection()[info.blockID];
//        
//        double h       = info.h[0];
//        double h2      = h*h;
//        Real m, mS, mN, mW, mE, mF, mB;
//        
//        printf("Loop through \n");
//        for(int iz=0; iz<B::sizeZ; ++iz)
//            for(int iy=0; iy<B::sizeY; ++iy)
//                for(int ix=0; ix<B::sizeX; ++ix)
//                {
//                    // get global index
//                    const int gix = ix + info.index[0] * B::sizeX;
//                    const int giy = iy + info.index[1] * B::sizeY;
//                    const int giz = iz + info.index[2] * B::sizeZ;
//                    
//                    const int idx = gix*7 + giy * 7 * blocksPerDimension * B::sizeX
//                    + giz * 7 * blocksPerDimension * B::sizeX * blocksPerDimension*B::sizeY;
//                    
//                    values[idx  ] =   4;
//                    values[idx+1] = - 1;
//                    values[idx+2] = - 1;
//                    values[idx+3] = - 1;
//                    values[idx+4] = - 1;
//                    values[idx+5] = - 1;
//                    values[idx+6] = - 1;
//                }
//        
//        printf("Loop through done, exiting matrix\n");
//        
//    }
    
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
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypre_grid, &rhs);
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypre_grid, &solution);
        
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
    
    
//    /* Set up Struct Vectors for rhs and solution */
//    inline void _setupVectors()
//    {
//        /* Create an empty vector object */
//        HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypre_grid, &rhs);
//        HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypre_grid, &solution);
//        
//        /* Indicate that the vector coefficients are ready to be set */
//        HYPRE_StructVectorInitialize(rhs);
//        HYPRE_StructVectorInitialize(solution);
//        
//        /* Set the vector coefficients */
//        const int nvalues      = B::sizeX * B::sizeY * B::sizeZ;
//        double * valuesRhs     = new double[nvalues];
//        double * valuesSol     = new double[nvalues];
//        
//        vector<BlockInfo> vInfo = mrag_grid->getBlocksInfo();
//        
//        printf("Hello my rank is %d \n", rank);
//        
//        
//        printf("Hello again my rank is %d \n", rank);
//        
//        //#pragma omp parallel for schedule(static)
//        for(int i = 0; i < blocksPerProcesor; i++)
//        {
//            int bID         = rank + (blocksPerProcesor - 1) * rank + i;
//            BlockInfo& info = vInfo[bID];
//            B& block        = mrag_grid->getBlockCollection()[info.blockID];
//            
//            int ilower[3]   = {info.index[0] * B::sizeX,
//                info.index[1] * B::sizeY,
//                info.index[2] * B::sizeZ };
//            
//            int iupper[3]   = { ilower[0] + B::sizeX - 1,
//                ilower[1] + B::sizeY - 1,
//                ilower[2] + B::sizeZ - 1};
//            
//            double h       = info.h[0];
//            double h2      = h*h;
//            
//            int idx = 0;
//            
//            for(int iz=0; iz<B::sizeZ; iz++)
//                for(int iy=0; iy<B::sizeY; iy++)
//                    for(int ix=0; ix<B::sizeX; ix++)
//                    {
//                        valuesRhs[idx] = h2 * block(ix,iy,iz).f;
//                        valuesSol[idx] =      block(ix,iy,iz).p;
//                        idx++;
//                    }
//            
//            HYPRE_StructVectorSetBoxValues(rhs, ilower, iupper, valuesRhs);
//            HYPRE_StructVectorSetBoxValues(solution, ilower, iupper, valuesSol);
//        }
//        
//        
//        /* This is a collective call finalizing the vector assembly.
//         The vectors are now ``ready to be used'' */
//        HYPRE_StructVectorAssemble(rhs);
//        HYPRE_StructVectorAssemble(solution);
//        
//        delete [] valuesRhs;
//        delete [] valuesSol;
//
//    }
//
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
    
    void _solveSystem()
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
    
//    void _getResultsOMP()
//    {
//        // get solution
//        int ilower[3]       = { 0, 0, 0 };
//        int iupper[3]       = { GridsizeX-1, GridsizeY-1, GridsizeZ-1 };
//        const int size3     = GridsizeX * GridsizeY * GridsizeZ;
//        double * values     = new double[size3];
//        HYPRE_StructVectorGetBoxValues(solution, ilower, iupper, values);
//        
//        vector<BlockInfo> vInfo = mrag_grid->getBlocksInfo();
//        
//#pragma omp parallel for schedule(static)
//        for(int i=0; i<vInfo.size(); i++)
//        {
//            const BlockInfo info = vInfo[i];
//            B& block = mrag_grid->getBlockCollection()[info.blockID];
//            
//            for(int iz=0; iz<B::sizeZ; iz++)
//                for(int iy=0; iy<B::sizeY; iy++)
//                    for(int ix=0; ix<B::sizeX; ix++)
//                    {
//                        const int gix = ix + info.index[0] * B::sizeX;
//                        const int giy = iy + info.index[1] * B::sizeY;
//                        const int giz = iz + info.index[2] * B::sizeZ;
//                        
//                        const int idx = gix + giy * blocksPerDimension * B::sizeX
//                        + giz * blocksPerDimension * B::sizeX * blocksPerDimension*B::sizeY;
//                        
//                        block(ix,iy,iz).p = values[idx];
//                    }
//        }
//        
//        delete [ ] values;
//    }
    
    
//    void _getResultsOMP()
//    {
//        vector<BlockInfo> vInfo = mrag_grid->getBlocksInfo();
//        
//        
//        if(rank == 0)
//        {
//            printf("OK\n");
//            //        for(int i = 0; i < blocksPerProcesor; i++)
//            for(int i = 0; i < vInfo.size(); i++)
//            {
//                int blockID     = i; //rank + (blocksPerProcesor - 1) * rank + i;
//                BlockInfo& info = vInfo[blockID];
//                B& block        = mrag_grid->getBlockCollection()[info.blockID];
//                
//                int ilower[3]   = {
//                    info.index[0] * B::sizeX,
//                    info.index[1] * B::sizeY,
//                    info.index[2] * B::sizeZ };
//                
//                int iupper[3]   = {
//                    ilower[0] + B::sizeX - 1,
//                    ilower[1] + B::sizeY - 1,
//                    ilower[2] + B::sizeZ - 1};
//                
//                int nvalues         = B::sizeX * B::sizeY * B::sizeZ;
//                double * values     = new double  [nvalues];
//                
//                HYPRE_StructVectorGetBoxValues(solution, ilower, iupper, values);
//                
//                
//                int idx             = 0;
//                
//                for(int iz=0; iz<B::sizeZ; iz++)
//                    for(int iy=0; iy<B::sizeY; iy++)
//                        for(int ix=0; ix<B::sizeX; ix++)
//                        {
//                            block(ix,iy,iz).p = values[idx];
//                            idx++;
//                        }
//                delete [] values;
//                
//            }
//        }
//    }
    
    
    
    void _getResultsOMP()
    {
        /* 1) Each processor fills results into blocks it owns, using own mrag_grid, since MPI doesn't share memory */
        vector<BlockInfo> vInfo = mrag_grid->getBlocksInfo();
        const int N = vInfo.size();
        
        for(int i = 0; i < blocksPerProcesor; i++)
        {
            int blockID     = rank + (blocksPerProcesor - 1) * rank + i;
            BlockInfo& info = vInfo[blockID];
            B& block        = mrag_grid->getBlockCollection()[info.blockID];
            
            int ilower[3]   = {
                info.index[0] * B::sizeX,
                info.index[1] * B::sizeY,
                info.index[2] * B::sizeZ };
            
            int iupper[3]   = {
                ilower[0] + B::sizeX - 1,
                ilower[1] + B::sizeY - 1,
                ilower[2] + B::sizeZ - 1};
            
            int nvalues         = B::sizeX * B::sizeY * B::sizeZ;
            double * values     = new double  [nvalues];
            
            HYPRE_StructVectorGetBoxValues(solution, ilower, iupper, values);
            
            int idx = 0;
            
            for(int iz=0; iz<B::sizeZ; iz++)
                for(int iy=0; iy<B::sizeY; iy++)
                    for(int ix=0; ix<B::sizeX; ix++)
                    {
                        block(ix,iy,iz).p = values[idx];
                        idx++;
                    }
            
            delete [] values;
            
        }
        
        /*2) Collect all the results into the mrag_grid owned by master rank */
        MPI_Request request[N];
        MPI_Status  status[N];
        
        if(nprocs > 1)
        {
            if(rank>0)
            {
                // - send
                for(int i = 0; i < blocksPerProcesor; i++){
                    int blockID     = rank + (blocksPerProcesor - 1) * rank + i;
                    BlockInfo& info = vInfo[blockID];
                    B& block        = mrag_grid->getBlockCollection()[info.blockID];
                    
                    int ierror = MPI_Isend(&block, sizeof(block), MPI_BYTE, 0, blockID, MPI_COMM_WORLD, &request[blockID]);
                }
            }
            else
            {
                // recieve
                for(int prcs=1; prcs<nprocs; prcs++)
                    for(int i = 0; i < blocksPerProcesor; i++)
                    {
                        int blockID     = prcs + (blocksPerProcesor - 1) * prcs + i;
                        BlockInfo& info = vInfo[blockID];
                        B& block        = mrag_grid->getBlockCollection()[info.blockID];
                        
                        int ierror = MPI_Recv(&block, sizeof(block), MPI_BYTE, prcs, blockID, MPI_COMM_WORLD, &status[blockID]);
                    }
            }
        }
    }



public:
    HelmholtzSolver3D_Hypre_MPI( ): bAlreadyAllocated(false)
    {  }
    
    ~HelmholtzSolver3D_Hypre_MPI()
    { }
    
    void setup_hypre(Grid<W,B>& input_grid, int rank, int nprocs, bool bVerbose=false, bool bCG=false, bool bRelaxation=false, std::vector<Real>* kappa = NULL, bool bMobility=false, std::vector<Real>* mobility = NULL)
    {
        mrag_grid           = &input_grid;
        this->rank          = rank;
        this->nprocs        = nprocs;
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
        
        this->GridsizeX         = blocksPerDimension * B::sizeX;
        this->GridsizeY         = blocksPerDimension * B::sizeY;
        this->GridsizeZ         = blocksPerDimension * B::sizeZ;
        this->nMRAGblocks       = blocksPerDimension * blocksPerDimension * blocksPerDimension;
        this->blocksPerProcesor = nMRAGblocks / nprocs;
        
        if((rank==0)&&(bVerbose)){
            printf("----------------------------------------------------------------------\n");
            printf("Hello from Helmholtz-Hypre solver \n");
            printf("Assume: uniform grid, which fits in the memory of single processor    \n");
            printf("Number of MRAG blocks is multiple of number of MPI processes.         \n");
            printf("Using %d MRAG blocks, %d MPI processes, %d blocksPerProcesor          \n", nMRAGblocks, nprocs, blocksPerProcesor);
            printf("----------------------------------------------------------------------\n");
        }
        
        assert( (nMRAGblocks % nprocs) == 0);
        
        if (bAlreadyAllocated)
            clean();
        
        _setup_grid();
        _setup_stencil();
        _setup_solver();
        
        bAlreadyAllocated = true;
        if ((bVerbose)&&(rank==0)) printf("Done with hypre setup!\n"); //exit(0);
    }
    
    
    void inline solve()
    {
        _setupMatrix();
        _setupVectors();
        _solveSystem();
        _getResultsOMP();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
//    void inline cleanUp()
//    {
//        HYPRE_StructGridDestroy(hypre_grid);
//        HYPRE_StructStencilDestroy(stencil);
//        HYPRE_StructMatrixDestroy(matrix);
//        HYPRE_StructVectorDestroy(rhs);
//        HYPRE_StructVectorDestroy(solution);
//        
//        if(bCG==0)
//            HYPRE_StructSMGDestroy(solver);
//        else
//        {
//            HYPRE_StructPCGDestroy(solver);
//            HYPRE_StructSMGDestroy(precond);
//        }
//        
//        bAlreadyAllocated = false;
//    }
    
    void inline clean()
    {
        HYPRE_StructGridDestroy(hypre_grid);
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
