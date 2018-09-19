#pragma once

#include <vtkPoints.h> 
#include <vtkCell.h>
#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

namespace MRAG
{
	
	template<typename TWavelets, typename TBlock, int nChannels=1, int iChannelStart=0>
	class IO_VTKNative {
	public:
		typedef MRAG::Grid<TWavelets, TBlock> GridType;
		typedef typename TBlock::ElementType ElementType;
		
		void Write(GridType & inputGrid, BoundaryInfo & bInfo,  string fileName );
	};
	
	
	template<typename TWavelets, typename TBlock, int nChannels, int iChannelStart>
	void
	IO_VTKNative< TWavelets, TBlock, nChannels, iChannelStart >::
	Write( GridType & inputGrid, BoundaryInfo & bInfo, string fileName_ )
	{
		// Open output file
		string filename  = fileName_ + ".vtu";
		// Calculate total number of points and cells
		vector<BlockInfo> vInfo = inputGrid.getBlocksInfo();
		unsigned int totalNumberOfPoints = 0;
		unsigned int verticesPerCell = 0;
		unsigned int totalNumberOfCells = 0;
		if( TBlock::sizeZ == 1 )
		{
			totalNumberOfPoints = vInfo.size() * (TBlock::sizeX+1) * (TBlock::sizeY+1);
			totalNumberOfCells = vInfo.size() * (TBlock::sizeX) * (TBlock::sizeY);
			verticesPerCell = 4;
		}
		
		// Preapare header
		vtkPoints * points = vtkPoints::New();
		points->SetNumberOfPoints(totalNumberOfPoints);
		float *pts = (float*) points->GetVoidPointer(0);
		for(int i=0; i<vInfo.size(); i++)
		{
			BlockInfo& info = vInfo[i];
			
			if( TBlock::sizeZ == 1 )
			{
				int iz = 0;
				int iPoint = 0;
				for(int iy=0; iy<TBlock::sizeY+1; iy++)
					for(int ix=0; ix<TBlock::sizeX+1; ix++)
					{
						float x[3];
						info.pos(x, ix, iy,iz);
						*pts++ = x[0] - TWavelets::CenteringOffset * info.h[0];
						*pts++ = x[1] - TWavelets::CenteringOffset * info.h[1];
						*pts++ = 0.0f;
						++iPoint;
					}
			}
		}
		
		vtkUnstructuredGrid * uGrid = vtkUnstructuredGrid::New();
		uGrid->SetPoints(points);
		uGrid->Allocate(totalNumberOfCells);
		vtkIdType verts[verticesPerCell];
		
		unsigned int counter = 0;
		for(int i=0; i<vInfo.size(); i++)
		{
			if( TBlock::sizeZ == 1 )
			{
				for(int iy=0; iy<TBlock::sizeY; iy++)
					for(int ix=0; ix<TBlock::sizeX; ix++)
					{
						for(int j = 0; j < verticesPerCell; j++)
						{
							int shift[3] = { j&1, (j>>1)&1, (j>>2)&1 };
							int ixx = ix + shift[0];
							int iyy = iy + shift[1];
							int izz = 0;
							int pointIndex = counter*(TBlock::sizeX+1)*(TBlock::sizeY+1)*(TBlock::sizeZ) + izz*(TBlock::sizeX+1)*(TBlock::sizeY+1) + iyy*(TBlock::sizeX+1) + ixx;
							verts[j] = pointIndex;
						}
						uGrid->InsertNextCell(VTK_PIXEL, verticesPerCell, verts);
					}
				counter += 1;
			}
		}
		
		for(int ichannel = 0; ichannel<nChannels; ++ichannel){
			//
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
			
			if( bIsCellCentered ){
				BlockLab<TBlock> lab;
				
				int steStart[3] ={0,0,0};
				int steEnd[3] = {2,2,1};
				
				lab.prepare(inputGrid.getBlockCollection(), bInfo ,steStart,steEnd);
				
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
					lab.load(info);
					int icount = 0;
					Real h = info.h[0];
					//using the blocklab, don't need this: \TBlock& block = inputGrid.getBlockCollection()[info.blockID];
					if( TBlock::sizeZ == 1 )
					{
						int iz = 0;
						//not used: unsigned int iz = 0;
						for(int iy=0; iy<TBlock::sizeY; iy++)
							for(int ix=0; ix<TBlock::sizeX; ix++)
							{
								const float rValue = ( lab(ix,iy,iz).giveMe(ichannel + iChannelStart, h) ) ;
								*vpts++ = rValue;
								++icount;
							}
						
					}
				}
				
				uGrid->GetCellData()->AddArray(fa);
				
			} else {
				BlockLab<TBlock> lab;
				
				int steStart[3] ={0,0,0};
				int steEnd[3] = {2,2,1};
				
				lab.prepare(inputGrid.getBlockCollection(), bInfo ,steStart,steEnd);
				
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
					lab.load(info);
					Real h = info.h[0];
					int icount = 0;
					//using the blocklab, don't need this: \TBlock& block = inputGrid.getBlockCollection()[info.blockID];
					if( TBlock::sizeZ == 1 )
					{
						int iz = 0;
						//not used: unsigned int iz = 0;
						for(int iy=0; iy<TBlock::sizeY+1; iy++)
							for(int ix=0; ix<TBlock::sizeX+1; ix++)
							{
								const float rValue = ( lab(ix,iy,iz).giveMe(ichannel + iChannelStart, h) ) ;
								*vpts++ = rValue;
								++icount;
							}
						
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
		
		points->Delete();
		uGrid->Delete();
		writer->Delete();
	}
}
