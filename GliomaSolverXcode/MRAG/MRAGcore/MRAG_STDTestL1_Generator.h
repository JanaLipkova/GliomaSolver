/*
 *  MRAG_STDTestL1_Generator.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/24/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "MRAG_STDTestL1.h"
#include "MRAG_STDTestL1_memtest.h"
#include "MRAG_STDTestL1_Grid.h"
#include "MRAG_STDTestL1_Refinement.h"
#include "MRAG_STDTestL1_Compression.h"
#include "MRAG_STDTestL1_BlockCollection.h"
#include "MRAG_STDTestL1_BoundaryInfo.h"

#include "MRAGWavelets_Interp2ndOrder.h"
#include "MRAGWavelets_Interp4thOrder.h"


namespace MRAG
{
class MRAG_STDTestL1_Generator
{
	static int argc;
	static char ** argv;
	
	template<typename Wavelets, int nX, int nY, int nZ>
	inline static void _test()
	{
	//	printf("Testing %d %d %d\n", nX, nY, nZ);
	//	MRAG_STDTestL1< Wavelets, Block<float, nX,nY,nZ> >::runTests();
	//	MRAG_STDTestL1_BlockCollection< Wavelets, Block<float, nX,nY,nZ> >::runTests();
	//	MRAG_STDTestL1_memtest< Wavelets, Block<float, nX,nY,nZ> >::runTests(2000,20);
		MRAG_STDTestL1_BoundaryInfo< Wavelets,  Block<float, nX,nY,nZ> >::runTests(argc, argv, false);
	}
	
	template<typename Wavelets, int n>
	class TEST
	{
	public:
		static void generate()
		{
			_test<Wavelets, n, n, n>();
			//_test<Wavelets, n, n, 1>();
			//_test<Wavelets, n, 1, 1>();
			
			
			TEST<Wavelets, n/2>::generate();
		}
	};
	
	template <typename Wavelets>
	class TEST< Wavelets, 8> { public:static void generate(){
		printf("End of the Testing\n");
	}};
	
	
public:
	static void run(int c, char ** v)
	{
		argc  = c;
		argv = v;
		TEST<Wavelets_Interp2ndOrder,32>::generate();
		TEST<Wavelets_Interp4thOrder,32>::generate();
	}
};
	
	int MRAG_STDTestL1_Generator::argc = 0;
	char** MRAG_STDTestL1_Generator::argv = NULL;
}