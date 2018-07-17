/*
 *  MRAGAddressable_CUDA.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/14/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <vector>
using namespace std;

namespace MRAG
{
	class Addressable_CUDA
	{
	public:
		virtual void * addressOf(const int blockID) const =0;
	};
	
}
