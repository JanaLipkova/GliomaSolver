/*
 *  MRAG_STDTestL1_Generator.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 7/24/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once


#include "MRAG_STDTestL2_Boundary.h"
#include "MRAG_STDTestL2_Compression.h"
#include "MRAG_STDTestL2_Refinement.h"
#include "MRAG_STDTestL2_IC.h"

#include "MRAGWavelets_Interp2ndOrder.h"
#include "MRAGWavelets_Interp4thOrder.h"


namespace MRAG
{
class MRAG_STDTestL2_Generator
{
	static int argc;
	static char ** argv;
	
	template<typename Wavelets, int nX, int nY, int nZ>
	inline static void _test()
	{
		printf("Testing %d %d %d\n", nX, nY, nZ);
		MRAG_STDTestL2_Boundary<MRAG::Wavelets_Interp2ndOrder, MRAG::Block<misc::memNum<float>, nX,nY,nZ> >::runTests(2,4,2,1);
		MRAG_STDTestL2_Compression< Wavelets, Block<float, nX,nY,nZ> >::runTests();
		MRAG_STDTestL2_Refinement< Wavelets, Block<float, nX,nY,nZ> >::runTests();
		MRAG_STDTestL2_IC< Wavelets,  Block<float, nX,nY,nZ> >::runTests();
	}
	
	template<typename Wavelets, int n>
	class TEST
	{
	public:
		static void generate()
		{
			//_test<Wavelets, n, n, n>();
			_test<Wavelets, n, n, 1>();
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
	
	int MRAG_STDTestL2_Generator::argc = 0;
	char** MRAG_STDTestL2_Generator::argv = NULL;
}