/*
 *  MPCFVtkIO.h
 *  MPCFTR
 *
 *  Created by Babak Hejazialhosseini  on 5/12/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <vtkPoints.h> 
#include <vtkCell.h>
#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

using namespace MRAG;

template<typename TWavelets, typename TBlock, int nChannels=1, int iChannelStart=0, typename TLab=BlockLab<TBlock> >
class IO_VTKNative3D
{
public:
	typedef MRAG::Grid<TWavelets, TBlock> GridType;
	typedef typename TBlock::ElementType ElementType;
	
	void Write(GridType & inputGrid, BoundaryInfo & bInfo,  string fileName );
};

template<typename TWavelets, typename TBlock, int nChannels, int iChannelStart, typename TLab>
void
IO_VTKNative3D< TWavelets, TBlock, nChannels, iChannelStart, TLab >::
Write( GridType & inputGrid, BoundaryInfo & bInfo, string fileName )
{
	string filename  = fileName + ".vtu";
	vector<BlockInfo> vInfo = inputGrid.getBlocksInfo();
	unsigned int totalNumberOfPoints = 0;
	unsigned int verticesPerCell = 0;
	unsigned int totalNumberOfCells = 0;
	
	totalNumberOfPoints = vInfo.size() * (TBlock::sizeX+1) * (TBlock::sizeY+1) * (TBlock::sizeZ+1);
	totalNumberOfCells = vInfo.size() * (TBlock::sizeX) * (TBlock::sizeY) * (TBlock::sizeZ);
	verticesPerCell = 8;
	
	vtkPoints * points = vtkPoints::New();
	points->SetNumberOfPoints(totalNumberOfPoints);
	float *pts = (float*) points->GetVoidPointer(0);
	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo& info = vInfo[i];
		
		int iPoint = 0;
		for(int iz=0; iz<TBlock::sizeZ+1; iz++)
			for(int iy=0; iy<TBlock::sizeY+1; iy++)
				for(int ix=0; ix<TBlock::sizeX+1; ix++)
				{
					float x[3];
					info.pos(x, ix, iy,iz);
					*pts++ = x[0] - TWavelets::CenteringOffset * info.h[0];
					*pts++ = x[1] - TWavelets::CenteringOffset * info.h[1];
					*pts++ = x[2] - TWavelets::CenteringOffset * info.h[2];
					++iPoint;
				}
	}
	
	vtkUnstructuredGrid * uGrid = vtkUnstructuredGrid::New();
	uGrid->SetPoints(points);
	points->Delete();
	uGrid->Allocate(totalNumberOfCells);
	vtkIdType verts[verticesPerCell];
	
	unsigned int counter = 0;
	for(int i=0; i<vInfo.size(); i++)
	{
		for(int iz=0; iz<TBlock::sizeZ; iz++)
			for(int iy=0; iy<TBlock::sizeY; iy++)
				for(int ix=0; ix<TBlock::sizeX; ix++)
				{
					for(int j = 0; j < verticesPerCell; j++)
					{
						int shift[3] = { j&1, (j>>1)&1, (j>>2)&1 };
						int ixx = ix + shift[0];
						int iyy = iy + shift[1];
						int izz = iz + shift[2];
						int pointIndex = counter*(TBlock::sizeX+1)*(TBlock::sizeY+1)*(TBlock::sizeZ+1) + izz*(TBlock::sizeX+1)*(TBlock::sizeY+1) + iyy*(TBlock::sizeX+1) + ixx;
						verts[j] = pointIndex;
					}
					uGrid->InsertNextCell(VTK_VOXEL, verticesPerCell, verts);
				}
		counter += 1;
	}
	
	for(int ichannel = 0; ichannel<nChannels; ++ichannel)
	{
		const bool bIsCellCentered = (TWavelets::CenteringOffset>0.0f);
		char channelName[256];
		sprintf(channelName,"channel%d",ichannel);
		vtkFloatArray * fa = vtkFloatArray::New();
		fa->SetNumberOfComponents(1);
		fa->SetName(channelName);
		if( bIsCellCentered )
			fa->SetNumberOfTuples(totalNumberOfCells);
		else
			fa->SetNumberOfTuples(totalNumberOfPoints);
		
		float * vpts = (float*)fa->GetVoidPointer(0);
		
		if( bIsCellCentered )
		{
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				TBlock& block = inputGrid.getBlockCollection()[info.blockID];
				
				int icount = 0;
				Real h = info.h[0];
				for(int iz=0; iz<TBlock::sizeZ; iz++)
					for(int iy=0; iy<TBlock::sizeY; iy++)
						for(int ix=0; ix<TBlock::sizeX; ix++)
						{
							const float rValue = ( block(ix,iy,iz).giveMe(ichannel + iChannelStart, h) ) ;
							*vpts++ = rValue;
							++icount;
						}
			}
			
			uGrid->GetCellData()->AddArray(fa);
			
		}
		else
		{
			TLab lab;
			
			int steStart[3] ={0,0,0};
			int steEnd[3] = {2,2,2};
			
			lab.prepare(inputGrid.getBlockCollection(), bInfo ,steStart, steEnd);
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				lab.load(info);
				Real h = info.h[0];
				int icount = 0;
				for(int iz=0; iz<TBlock::sizeZ+1; iz++)
					for(int iy=0; iy<TBlock::sizeY+1; iy++)
						for(int ix=0; ix<TBlock::sizeX+1; ix++)
						{
							const float rValue = ( lab(ix,iy,iz).giveMe(ichannel + iChannelStart, h) ) ;
							*vpts++ = rValue;
							++icount;
						}
			}
			uGrid->GetPointData()->AddArray(fa);
		}
		fa->Delete();
	}
	
	vtkXMLUnstructuredGridWriter * writer = vtkXMLUnstructuredGridWriter::New();
	writer->SetFileName(filename.c_str());
	writer->SetInput(uGrid);
	writer->Write();
	
	uGrid->Delete();
	writer->Delete();
}