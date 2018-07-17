/*
 *  MRAGVisualTypes.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 9/19/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once
namespace MRAG
{
	namespace Visual
	{
		
		struct RGBA
		{
			float r,g,b,a;
			
			RGBA():r(0), g(0), b(0), a(0){}
			RGBA(float x,float y,float z,float w):r(x), g(y), b(z), a(w){}
			
			RGBA(const RGBA& r):r(r.r), g(r.g), b(r.b), a(r.a){}
		};
		
		template <typename T> 
		RGBA convertToRGBA(const T& p)
		{
			return RGBA();
		}
		
	}
}
