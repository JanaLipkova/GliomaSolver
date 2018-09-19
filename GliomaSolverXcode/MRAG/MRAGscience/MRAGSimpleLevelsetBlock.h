/*
 *  MRAGSimpleLevelsetBlock.h
 *  MRAG
 *
 *  Created by Michael Bergdorf on 9/9/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGrid.h"


namespace MRAG {
   
   template <typename DataType, int cNarrowBand =4 , int cSizeX=16, int cSizeY=1, int cSizeZ=1>
   struct SimpleLevelsetBlock {
   typedef DataType ElementType;
   
	DataType data[cSizeZ][cSizeY][cSizeX];
   
	SimpleLevelsetBlock(DataType v = DataType()):h(0), bHIsSet(false)
	{
      if(cDoDebugChecks) assert(cNarrowBand>1);
		for(int iz=0; iz<cSizeZ; iz++)
			for(int iy=0; iy<cSizeY; iy++)
				for(int ix=0; ix<cSizeX; ix++)
					data[iz][iy][ix] = v;
	}
	
	inline DataType& operator()(int ix, int iy=0, int iz=0)
	{
		if (cDoDebugChecks)
		{
			assert(ix>=0 && ix<cSizeX);
			assert(iy>=0 && iy<cSizeY);
			assert(iz>=0 && iz<cSizeZ);
			assert((ElementType *)&data[iz][iy][ix] == (ElementType*)data + cSizeX*cSizeY*iz + cSizeX*iy + ix);
		}
		
		return data[iz][iy][ix];
	}
	
	inline DataType& operator[](int i)
	{
		if (cDoDebugChecks)
		{
			assert(i>=0);
			assert(i<cSizeX*cSizeY*cSizeZ);
		}
		
		return *((ElementType *)data+i);
	}
	
	inline void wtdata(Real& coeff, int code, int ix, int iy=0, int iz=0)
	{
      if(cDoDebugChecks){
        // assert(bHIsSet);
      }
	  
	  if (!bHIsSet) return;
      Real sLevelSetValue = data[iz][iy][ix].levelset()/h;
      const Real sClippingValue = _eta(sLevelSetValue,cNarrowBand);
      coeff *= sClippingValue;
      return;
	}
   
   void setH(Real _h){
      this->h = _h;
      this->bHIsSet = true;
   }
	
   
   template<typename T, int i> 
   static inline Real levelset_projector_impl(const T&t) 
   {
      return (Real)t.levelset();
   }
   
   /* Level set specific */
   Real h;
   bool bHIsSet;
      
	static const int sizeX = cSizeX;
	static const int sizeY = cSizeY;
	static const int sizeZ = cSizeZ;
   
	static const bool shouldProcessDirectionX = cSizeX>1; 
	static const bool shouldProcessDirectionY=  cSizeY>1; 
	static const bool shouldProcessDirectionZ = cSizeZ>1; 
	
	bool shouldBeRefined() { return false; }
	bool shouldBeCompressed() { return false; }
	inline static int getMemSize()  { return cSizeX*cSizeY*cSizeZ*sizeof(DataType); }
   
private:
   inline const Real _eta(Real phi, const int maxPhi) const {
      //
      const Real gamma = (Real)maxPhi;
      const Real beta  = gamma-1.0;
      const Real phimapped  = fabs(phi); 
      Real eta;
      if(phimapped<gamma){
         if(phimapped<(gamma-1.0)){
            eta = 1.0;
         } else {
            eta = pow((double)phimapped-gamma,2.0)*(2.0*phimapped+gamma-3.0*beta);
         }
      } else {
         eta = 0.0;
      }
      return eta;
   }
};

}

