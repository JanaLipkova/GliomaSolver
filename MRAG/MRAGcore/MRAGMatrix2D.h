/*
 *  Matrix3D.h
 *  ReactionDiffusion
 *
 *  Created by Diego Rossinelli on 10/19/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#pragma once
#include <iostream>
using namespace std;

namespace MRAG
{

template <class DataType>  class Matrix2D
{
private:
	DataType * m_pData;
	unsigned int m_vSize[2];
	unsigned int m_nElements;

public:
	Matrix2D(unsigned int nSizeX, unsigned int nSizeY);
	~Matrix2D() { delete [] m_pData; }
	
	DataType & Access(unsigned int ix, unsigned int iy) const;
	DataType & LinAccess(unsigned int i) const;
	
	unsigned int * getSize();
	unsigned int getNumberOfElements() const;
	
	Matrix2D& operator= (DataType val);
};

template <class DataType>  Matrix2D<DataType>::Matrix2D(unsigned int nSizeX, unsigned int nSizeY):
	m_pData(NULL),
	m_nElements(0)
{
	m_vSize[0] = nSizeX;
	m_vSize[1] = nSizeY;
	
	m_nElements = nSizeX*nSizeY;
	
	m_pData = new DataType[m_nElements];
}

template <class DataType> inline DataType & Matrix2D<DataType>::Access(unsigned int ix, unsigned int iy) const
{
	ix = (ix + m_vSize[0]) % m_vSize[0];
	iy = (iy + m_vSize[1]) % m_vSize[1];
	
#ifdef _DEBUG
	if (!(ix<m_vSize[0]) || !(iy<m_vSize[1]) ) printf("excpetion %d %d \n", ix,iy);
	assert(ix<m_vSize[0]);
	assert(iy<m_vSize[1]);
#endif
	
	return m_pData[ iy*m_vSize[0] + ix];
}

template <class DataType> inline DataType & Matrix2D<DataType>::LinAccess(unsigned int i) const
{
#ifdef _DEBUG
	assert(i<m_nElements);
#endif
	return m_pData[i];
}

template <class DataType> inline unsigned int  Matrix2D<DataType>::getNumberOfElements() const
{
	return m_nElements;
}

template <class DataType> inline unsigned int *  Matrix2D<DataType>::getSize() 
{
	return m_vSize;
}

template <class DataType> inline Matrix2D<DataType>&  Matrix2D<DataType>::operator= (DataType val)
{
	for(unsigned int i=0;i<m_nElements;i++) m_pData[i] = val;
	return *this;
}

}

