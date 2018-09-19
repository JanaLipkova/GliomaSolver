/*
 *  MRAG_STDTestL1_memtest.h
 *  MRAG
 *
 *  Created by Manfred on 7/18/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _MRAG_MEMTEST_
#define _MRAG_MEMTEST_

#include "MRAG_STDTestL1.h"
#include <iostream>
#include <string>
#include <sstream>
#include <sys/types.h>
#include <unistd.h>

namespace misc{

//little helper function: get memory usage by a call to ps -x -p <mypid> -o rss
//returns usage in MB
float getmemusage() {
    
	//TODO: make a system check (ifdef WIN_32)
    // or: during runtime
    float mbused;

	FILE *fp;    //output from popen
	const int buf_length=80;
	char psres[buf_length];

	std::string cur_string;
	std::vector<std::string> psretlist;
	std::ostringstream pscallstrm;
	std::string pscall;
	
	int kbused;
    pid_t mypid; //my process id
    
	mypid=getpid();
    
	std::cout << "Calling Memtest ... for process: " << mypid << std::endl ;
	pscallstrm << "ps -x -p " << mypid << " -o rss"; //returns size in kb (1024 byte)
	pscall=pscallstrm.str(); 
	
	fp=popen(pscall.c_str(),"r");
	if (fp !=NULL )
	{
		while(fgets(psres,buf_length,fp) != NULL)
		{
			cur_string=psres;
			if (cur_string [cur_string.size () - 1] != '\n')
			{
				std::cout <<std::endl << "could not interpretate result of ps-call. check buf_length" <<std::endl;
				break;
			}
			assert (cur_string [cur_string.size () - 1] == '\n');
			psretlist.push_back (cur_string.substr (0, cur_string.size () - 1));
		}
		pclose(fp);
		if ( psretlist.size()>0 && (psretlist[0].find("RSS", 0) != std::string::npos))
		{
			kbused=atoi(psretlist[1].c_str());
			mbused=float(kbused)/float(1024.0);
		}
		else
		{ 
			std::cout << "could interpretate result of ps-call: does ps exist on your system?" <<std::endl;
			mbused=-1.0;
		}
	}
	else
	{
		std::cout <<std::endl << "could not open pipe to process: does popen() exist on your system?" <<std::endl;
		mbused=-1.0;
	}
	
	//clear:
	cur_string.resize(0);
	psretlist.resize(0);
	pscallstrm.clear();
	pscall.resize(0);
	
	
	return mbused;

}

}

namespace MRAG {

template <typename Wavelets, typename Block>
class MRAG_STDTestL1_memtest: public MRAG_STDTestL1<Wavelets, Block>
{

protected:

	bool run(int nallocs, int nBlocks)
	{
	  //1.) prepare system calls
	  //2.) allocate
	  //3.) calculate memory usage estimate
	  //4.) get real memory usage (through system call)
	  //5.) compare
	  //6.) deallocate, goto 2.)
	  


	  bool res=true;
	  
	  const float tol=0.1;
	  
  	  
	  float mbself;
	  float mbparts[3];
	  float mbused(0.0),mbold(0.0),mbshould(0.0);
	  float* mbincs=new float[nallocs];
	  
    
	  //own usage: integers+floats (neglecting usage in getmemusage)
	  mbself=((6+2)*4.0+(nallocs+3+4)*4.0)/float(1024.0*1024.0);
	  
	  for (int i=0;i<nallocs;++i)
	  {
          
	      mbold=mbused;
		  mbused=0;
		  mbparts[0]=0.0;
		  mbparts[1]=0.0;
		  mbparts[2]=0.0;
		  mbshould=0;
		  
		  try
			{
				MRAG::Grid<Wavelets, Block> * grid = new MRAG::Grid<Wavelets, Block>(nBlocks,nBlocks,1,false);
       			mbparts[0]=grid->getMemorySize();
				mbshould+=mbparts[0];
				
				//add up: Collection and Boundarinfo
				const BoundaryInfo& bInfo=grid->getBoundaryInfo();
				mbparts[1]=bInfo.getMemorySize();
		    	const BlockCollection<Block>& bCollection = grid->getBlockCollection();
  			    mbparts[2]=bCollection.getMemorySize();
				mbshould+=mbparts[1]+mbparts[2];
				
				
				mbused=misc::getmemusage();
				
				if (mbused > 0.0 )
				  {
					 
					 mbshould+=mbself;
					 std::cout <<std::endl<< " -------------------------MEMORY USAGE REPORT: ------------------------"<<std::endl;
					 std::cout << "Real Usage: [MB]:" << mbused <<std::endl;
					 std::cout << "should be [MB]:" << mbshould <<std::endl;
					 std::cout << "parts: Self:" << mbself <<std::endl;
					 std::cout << "Grid:" << mbparts[0] <<std::endl;
					 std::cout << "Boundary Information:" << mbparts[1] <<std::endl;
					 std::cout << "Collection:" << mbparts[2] <<std::endl; 
					 std::cout << "relative mismatch: " << fabs(mbused/mbshould) << std::endl;
					 
				  }
				else
					{
					throw "could not determine Memory usage from System. [getmemusage returned -1.0] " ;
					break;
					}
				
     			delete grid;
				mbincs[i]=mbused-mbold;
			}
			catch(...)
			{
				return false;
			}
			
			if ((fabs(mbused/mbshould-1.0) <  tol)&&res)
				{
					res= true;
				}
			else
				{
					res=false;
				}
		}
		
		
		std::cout <<std::endl<< " ------------MEMORY GROWTH REPORT:  Increments -------------"<<std::endl;
		std::cout << "[ ";
		for (int i=0;i<nallocs;++i)
		{
			std::cout << mbincs[i] << " " ;
		}
        std::cout << "]" <<std::endl;
		
		delete mbincs;
		
		return res;

	}


public:
	static void runTests(int nallocs=1, int nBlocks=20)
	{
		MRAG_STDTestL1_memtest<Wavelets,Block> test;
		test.run(nallocs,nBlocks);
	}

};

}



#endif