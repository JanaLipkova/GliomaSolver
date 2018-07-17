/*
 *  MRAGMacros.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 9/3/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once
#ifndef _MRAGHelperMacros_
#define _MRAGHelperMacros_
#include "MRAGCommon.h"
namespace MRAG
{
    /**
     * Creates a Projector class that provides a static function Real Project(const T&).
     * Macro is to be used at the top level where you would also define a class.
     * To be used for BlockFWT where result of Project is used to guide refinements/compressions.
     * @param name      Name used for the created class
     * @param f         Global function f of the form template <typename T, int i> inline Real f(const T&t).
     *                  Provides the result returned by static function Real Project(const T&).
     *                  With i one can number a set of projectors.
     * @see MRAG::BlockFWT
     */
	#define make_projector(name, f) \
	class name \
	{ \
	public:\
		template <typename T, int i>\
		inline static Real Project(const T& t)\
		{\
			return f<T, i>(t);\
		}\
	};
		
	 template <typename T, int i> inline Real dummy_projector_impl(const T&t)
	{
		return *(Real *)&t;
	}
	
#ifndef _DUMMY_PROJECTOR
#define _DUMMY_PROJECTOR
	make_projector(dummy_projector, dummy_projector_impl);
#endif
}
#endif //MRAGHelperMacros
