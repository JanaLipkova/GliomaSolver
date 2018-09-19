/*
 *  MRAG_STDTestL2_IC.h
 *  MRAG
 *
 *  Created by basil bayati on 7/30/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#pragma once

#include "MRAG_STDTestL2.h"
#include <fstream>
#include <sstream>

// This class contains methods for the initial condition.  
// All methods must return a double and take 4 doubles as
// parameters
namespace MRAG
{
class ICClass
{
	public:
	
	double sinCosineFunction(double x, double y, double z, double t)
	{
		return sin(2.0*M_PI*x)*cos(2.0*M_PI*y);
	}
	
	double diskFunction(double x, double y, double z, double t)
	{
		double radius = 0.2;
		double s = 2.0;
		double b = 160.0 / sqrt(6.0);
		double a = sqrt(2.0) - 1.0;
		double alpha = sqrt( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) );
		alpha -= radius;
		return ( 1.0 / (  pow(1.0 + a*exp(b*alpha), s)  ) );
	}
	double noisyDiskFunction(double x, double y, double z, double t)
	{
		double radius = 0.2;
		double s = 2.0;
		double b = 160.0 / sqrt(6.0);
		double a = sqrt(2.0) - 1.0;
		double alpha = sqrt( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) );
		alpha -= radius;
		
		double value = ( 1.0 / (  pow(1.0 + a*exp(b*alpha), s)  ) );
		if ( value < 1.0e-4 )
			return 0.0;
		
		return std::max(0.0, ( value + (drand48() - 0.5)/10.0));
	}
};
}

namespace MRAG
{
	template <typename Wavelets, typename Block>
	class MRAG_STDTestL2_IC : public MRAG_STDTestL2< Wavelets, Block > 
	{
		private:	
			typedef double (MRAG::ICClass::*InitialConditionPointer)(double, double, double, double); 
			
			int frameNumber;
			double currentMinDetail;
			double currentMaxDetail;
			vector< double > detailCoefficients;
			double epsilon;
			
			// for writing out the convergence data
			vector< double > Convergence_L2_Error_Vector;
			vector< double > Convergence_NumPoint_Vector;
			
			// The _write* methods write a tab delimited ASCII file.
			// Mathematica:
			// data=Import["/Users/basil/Desktop/MRAG/build/Debug/output0.txt","TSV"];
			// ListPointPlot3D[data,PlotRange->{-1.3,1.3},ViewPoint->{1,1.2,1.2},ImageSize->{500,500}]
			// 
			void _write2DField()
			{
				FILE* myfile;
	
				stringstream outP;
				string fileBase = "output";
				outP << frameNumber;
				string frameS = outP.str();
				string fileNameEnd = ".txt";
				string fileName;
				
				fileName = fileBase + frameS;
				fileName = fileName + fileNameEnd;
				
				cout << "---------------writing file: " << fileName << " ---------------" << endl;
				
				myfile = fopen(fileName.c_str(), "w");
				
				// write the data
				vector<BlockInfo> vInfo = this->grid.getBlocksInfo();
				
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
					Block& block = this->grid.getBlockCollection()[info.blockID];
					
					const double h[2]= {pow(2.,-info.level), pow(2.,-info.level)};
					const double start[2] = {info.index[0]*h[0],info.index[1]*h[1]};
					const int n[2] = {Block::sizeX+1, Block::sizeY+1};
					const double d[2] = {h[0]/(n[0]-1),h[1]/(n[1]-1)};
					const double cOffset = Wavelets::CenteringOffset;
					
					for(int iy=0; iy<Block::sizeY; iy++)
						for(int ix=0; ix<Block::sizeX; ix++)
						{
							const double x[2] = {start[0]+(ix)*d[0], start[1]+(iy)*d[1]};
							fprintf(myfile, "%lf\t%lf\t%lf\n", (x[0]+cOffset), (x[1]+cOffset), block(ix,iy) );
						}
				}
				
				++frameNumber;
				fclose(myfile);
			}
			
			void _writeVectors()
			{
				FILE* myfile;
				string fileName = "convergence.txt";
							
				cout << "---------------writing file: " << fileName << " ---------------" << endl;
				
				myfile = fopen(fileName.c_str(), "w");
				
				// write the data
				for(int i = 0; i < Convergence_L2_Error_Vector.size(); ++i)
				{
					fprintf(myfile, "%lf\t%lf\n",  Convergence_NumPoint_Vector[i], Convergence_L2_Error_Vector[i] );
				}
			
				fclose(myfile);
			}
			
			void _ic( InitialConditionPointer icPointer2Function )
			{
				MRAG::ICClass initialConditionObj;
				
				// put the initial condition on the grid
				vector<BlockInfo> vInfo = this->grid.getBlocksInfo();
				
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
					Block& block = this->grid.getBlockCollection()[info.blockID];
					
					const double h[2]= {pow(2.,-info.level), pow(2.,-info.level)};
					const double start[2] = {info.index[0]*h[0],info.index[1]*h[1]};
					const int n[2] = {Block::sizeX+1, Block::sizeY+1};
					const double d[2] = {h[0]/(n[0]-1),h[1]/(n[1]-1)};
					const double cOffset = Wavelets::CenteringOffset;
					
					for(int iy=0; iy<Block::sizeY; iy++)
						for(int ix=0; ix<Block::sizeX; ix++)
						{
							const double x[2] = {start[0]+(ix)*d[0], start[1]+(iy)*d[1]};
							block(ix,iy) = (initialConditionObj.*icPointer2Function)( (x[0]+cOffset), (x[1]+cOffset), 0.0, 0.0);
						}
				}
			}	
			
			void _fwt()
			{
				BlockFWT<Wavelets, Block> fwt;
				fwt.prepare(this->grid.getBlockCollection(), this->grid.getBoundaryInfo());
				vector<MRAG::BlockInfo> vInfo = this->grid.getBlocksInfo();
				detailCoefficients.clear();
				
				double maxDetail = 0;
				double minDetail = (double)HUGE_VAL;
				for(int i=0; i<vInfo.size(); i++)
				{
					fwt.template fwt<0>(vInfo[i]);
					
					const FWTReport<>& report= fwt.getReport();
					maxDetail = std::max(maxDetail, report.getCoeffMaxMag<0>(1,1));
					minDetail = std::min(minDetail, report.getCoeffMinMag<0>(1,1));
					
					detailCoefficients.push_back(report.getCoeffMaxMag<0>(1,1));
					cout << "i in fwt: " << i << endl;
				}
				cout << "max detail: " << maxDetail << endl;
				cout << "min detail: " << minDetail << endl;
				
				this->currentMinDetail = minDetail;
				this->currentMaxDetail = maxDetail;
			}
		
			void _refinements() 
			{
				set<int> shouldBeRefined;
				vector<MRAG::BlockInfo> vInfo = this->grid.getBlocksInfo();
				
				for(int i=0; i<vInfo.size(); i++)
				{
					cout << "i in _refinements: " << i << endl;
					
					double currentDetail = detailCoefficients[i] / currentMaxDetail;
					cout << "dt c: " << detailCoefficients[i] << endl;
					cout << "currentDetail: " << currentDetail << endl << endl;
					if (currentDetail > epsilon) 
						shouldBeRefined.insert(vInfo[i].blockID);
				}
				
				this->grid.refine(shouldBeRefined);
			}
			
			void _init()
			{
				this->frameNumber		= 0;
				this->currentMinDetail	= (double)HUGE_VAL;
				this->currentMaxDetail	= 0.0;
				Convergence_L2_Error_Vector.clear();
				Convergence_NumPoint_Vector.clear();
			}
			
			void _error( InitialConditionPointer icPointer2Function, double & L2_error, double & numberOfElements  )
			{
				MRAG::ICClass initialConditionObj;
				
				// put the initial condition on the grid
				vector<BlockInfo> vInfo = this->grid.getBlocksInfo();
				L2_error			= 0.0;
				numberOfElements	= 0.0;
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
					Block& block = this->grid.getBlockCollection()[info.blockID];
					
					const double h[2]= {pow(2.,-info.level), pow(2.,-info.level)};
					const double start[2] = {info.index[0]*h[0],info.index[1]*h[1]};
					const int n[2] = {Block::sizeX+1, Block::sizeY+1};
					const double d[2] = {h[0]/(n[0]-1),h[1]/(n[1]-1)};
					const double cOffset = Wavelets::CenteringOffset;
					
					for(int iy=0; iy<Block::sizeY; iy++)
						for(int ix=0; ix<Block::sizeX; ix++)
						{
							numberOfElements += 1.0;
							const double x[2] = {start[0]+(ix)*d[0], start[1]+(iy)*d[1]};
							L2_error += pow( (block(ix,iy) - (initialConditionObj.*icPointer2Function)( (x[0]+cOffset), (x[1]+cOffset), 0.0, 0.0)), 2.0 ); 
						}
				}
				L2_error /= numberOfElements;
				L2_error = sqrt(L2_error);
			}
			
			bool run( double epsilon, int numberOfLevels )
			{
				try
				{
					this->epsilon = epsilon;
					cout << "Start Test with epsilon: " << epsilon << ", numberOfLevels: " << numberOfLevels << endl;
										
					// set a pointer to the initial condition function contained in the class "ICClass"
					InitialConditionPointer icPointer2Function = &MRAG::ICClass::diskFunction;
					
					_init();
					_ic( icPointer2Function );
					_write2DField();
					
					double l2Error		= 0.0;
					double numberPoints = 0.0;
						
					for (int i = 0; i < numberOfLevels; ++i )
					{
						_fwt();										cout << "_fwt			i: " << i << endl;
						_refinements();								cout << "_refinements	i: " << i << endl;
						
						_error( icPointer2Function, l2Error, numberPoints );
						Convergence_L2_Error_Vector.push_back( l2Error );
						Convergence_NumPoint_Vector.push_back( numberPoints );
						
						_ic( icPointer2Function );					cout << "_ic			i: " << i << endl;
						_write2DField();							cout << "_write2DField			i: " << i << endl;						
					}
					
					// write out the data
					_writeVectors();
									
					cout << "Start Test with epsilon: " << epsilon << ", numberOfLevels: " << numberOfLevels << endl;
				}
				catch(...)
				{
					return false;
				}
				
				return true;
			}
			
		public:
					
			static void runTests()
			{
				MRAG_STDTestL2_IC<Wavelets, Block> test;
				assert( test.run(1.e-6, 6) == true );
			}
		};
}

