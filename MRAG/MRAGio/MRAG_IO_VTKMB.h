#pragma once

#include <string>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkZLibDataCompressor.h>

#include "vtkCompositeDataIterator.h" 
#include "vtkMultiBlockDataSet.h" 
#include "vtkCompositeDataPipeline.h" 

using namespace MRAG;
using namespace std;

template<typename TWavelets, typename TBlock, int nChannels=1,  int iChannelStart=0, typename TLab=BlockLab<TBlock> >
class IO_VTKNative_MB {
public:
	typedef MRAG::Grid<TWavelets, TBlock> GridType;
	typedef typename TBlock::ElementType ElementType;

	void Write(GridType & inputGrid, BoundaryInfo & bInfo, string fileName)
	{
		char buffer[255];
		sprintf(buffer, "%s.vtm", fileName.c_str());
		ofstream fwrapper(buffer);

		// Write vtm header
		fwrapper << "<?xml version=\"1.0\"?>" << endl
				 <<	"<VTKFile type=\"vtkMultiBlockDataSet\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl
				 << "	<vtkMultiBlockDataSet>" << endl;
		
		// Write individual image files
		vector<BlockInfo> vInfo = inputGrid.getBlocksInfo();
		
		TLab lab;
		int steStart[3] ={0,0,0};
		int steEnd[3] = {2,2,2};
		lab.prepare(inputGrid.getBlockCollection(), bInfo , steStart, steEnd);

		for(int i=0; i<vInfo.size(); i++)
		{
			BlockInfo& info = vInfo[i];
			VTI_Dump(info, lab, fileName.c_str(), i, inputGrid.getCurrentMaxLevel());	
			fwrapper << "		<DataSet dataset=\"" << i << "\" group=\"0\" file=\"" << fileName << "_" << i << ".vti\"" << "/>" << endl;
		}
		
		// Close wrapper
		fwrapper << "	</vtkMultiBlockDataSet>" << endl
				 << "</VTKFile>" << endl;
		fwrapper.close();
	}
	
	void VTI_Dump(BlockInfo& info, TLab& lab, const char* fileroot, int bidx, int Lmax)
	{	
		// Open output file
		char filename[255];
		sprintf(filename, "%s_%d.vti", fileroot, bidx);
		
		// Calculate total number of points and cells
		unsigned int totalNumberOfPoints = (TBlock::sizeX+1) * (TBlock::sizeY+1) * (TBlock::sizeZ+1);
		
		// Get origin
		float x[3];
		info.pos(x, 0,0,0);

		// Create an image data
		vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
		imageData->SetExtent(0,TBlock::sizeX, 0,TBlock::sizeY, 0,TBlock::sizeZ);
		imageData->SetDimensions(TBlock::sizeX+1, TBlock::sizeY+1, TBlock::sizeZ+1);
		imageData->SetNumberOfScalarComponents(nChannels);
		imageData->SetScalarTypeToFloat();
		const float h = info.h[0];
		imageData->SetSpacing(h,h,h);
		
		Real p[3];
		info.pos(p, 0,0,0);
		imageData->SetOrigin(p[0], p[1], p[2]);
		
		lab.load(info);
		for(int ichannel = 0; ichannel<nChannels; ++ichannel)
		{
			char channelName[256];
			sprintf(channelName,"channel%d",ichannel);
			vtkFloatArray * fa = vtkFloatArray::New();
			fa->SetNumberOfComponents(1);
			fa->SetName(channelName);
			fa->SetNumberOfTuples(totalNumberOfPoints);
			float * vpts = (float*)fa->GetVoidPointer(0);
			
			for(int bz=0; bz<TBlock::sizeZ+1; bz++) 
				for(int by=0; by<TBlock::sizeY+1; by++) 
					for(int bx=0; bx<TBlock::sizeX+1; bx++)
					{					
						const float rValue = lab(bx,by,bz).giveMe(ichannel+iChannelStart, h); 
						*vpts++ = rValue;
					} 
			imageData->GetPointData()->AddArray(fa);
			fa->Delete();		
		}	
		
		// Write to file	
		vtkXMLImageDataWriter *imageWriter = vtkXMLImageDataWriter::New();
		imageWriter->SetFileName(filename);
		imageWriter->SetNumberOfPieces(1);
		imageWriter->SetDataModeToBinary();
		vtkZLibDataCompressor* compressor = vtkZLibDataCompressor::New();
		imageWriter->SetCompressor(compressor);
		imageWriter->SetInput(imageData);
		imageWriter->Write();
		imageWriter->Delete();
	}
	
	
};