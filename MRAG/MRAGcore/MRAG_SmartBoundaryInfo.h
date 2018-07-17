/*
 *  MRAG_SmartBoundaryInfo.h
 *  MRAG
 *
 *  Created by Manfred Quack on 9/29/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 *  initial Version: 0.1 (still todo: memory leak check)
 */


#include "MRAGCommon.h"
#include "MRAGrid.h"
#include "MRAGBoundaryBlockInfo.h"

namespace MRAG {


//SmartBoundaryInfo class:
// Derived from BoundaryInfo (in: "MRAGBoundaryBlockInfo.h")
// Purpose: If you require a larger stencil than internally provided by MRAG, or to change the stencil during runtime, 
//         you can use this class.
// It detects whether it can use the internal stencil, or whether it should create a new boundary-set with _createBoundaryInfo.
// --> use init(Grid<Wavelets,Block>& grid, const int steStart[3], const int steEnd[3]) to initialize the smartBoundaryInfos the first time.
//     you can also use init(..) later again to change the stencil during runtime.
// --> use update to update the SmartBoundaryInfo (this will internally either call getBoundaryInfo(), or createBoundaryInfo() )
// --> removal of stuff is handled (hopefully correctly)

class SmartBoundaryInfo: public BoundaryInfo
	{
	 enum BoundaryInfoStatus
	 {
		Constructed,
		initialized_linking,
		initialized_owning,
		updating_linking,
		updating_owning,
	 };
	 
	 
	 
	 protected:
	 static const int _verb=1;
	 BoundaryInfoStatus bStatus;
	 int _smartSteStart[3],_smartSteEnd[3]; //I need this for the update-function. (integers)
	 bool bOwnerOfBoundaryInfo;
	 BoundaryInfo* owned_bInfo;
	 BoundaryInfo* linked_bInfo;
	  
	 void _shallow_clear()
	 {
		boundaryInfoOfBlock.clear();
		for (int i=0; i<3; ++i)
			{
				stencil_start[i]=0;
				stencil_end[i]=0;
			}
	 }
	 
	void _updateBase( BoundaryInfo*bi ) 
	{
	   assert(bi!=NULL);
	   //update baseclass
		boundaryInfoOfBlock=bi->boundaryInfoOfBlock;
		for (int i=0; i<3; ++i)
			{
				stencil_start[i]=bi->stencil_start[i];
				stencil_end[i]=bi->stencil_end[i];
			}
	}
	
		
	 //update based on grid.createBoundaryInfo.
	 template <typename Wavelets, typename Block>
	 void _updateOwning(Grid<Wavelets,Block>& grid)// update when owning (createBoundaryInfo)
		{
		assert(bStatus==updating_linking || bStatus==initialized_owning);
		assert(owned_bInfo!=NULL && linked_bInfo ==NULL);
		//delete old data create new one.
		BoundaryInfo* tmpBnd=owned_bInfo;
		delete tmpBnd;
		if (_verb>0) std::cout << "Update-Owning: Calling createBoundaryInfo..." << std::endl;
		owned_bInfo=grid.createBoundaryInfo(_smartSteStart,_smartSteEnd);
		_updateBase(owned_bInfo);
 		}
	
 	 //update based on grid.getBoundaryInfo:
	 template <typename Wavelets, typename Block>
	 void _updateLinking(Grid<Wavelets,Block>& grid)  // update when using MRAG-Boundaries (getBoundaryInfo)
		{
		assert(bStatus==updating_linking || bStatus==initialized_linking);
		assert(owned_bInfo==NULL && linked_bInfo !=NULL);
		if (_verb>0) std::cout << "Update-Internal: Calling getBoundaryInfo..." <<std::endl;
		//just updated pointer to boundaryInfo.
		linked_bInfo=&grid.getBoundaryInfo();
		_updateBase(linked_bInfo);		
		}
	 

	 
	 public:
	 //SmartBoundaryInfo(BoundaryInfo* bInfo, bCreated=true)
	 //SmartBoundaryInfo(BoundaryInfo& bInfo)

   
	 
	 //Constructor:
	 SmartBoundaryInfo():
	 bOwnerOfBoundaryInfo(false),bStatus(Constructed),owned_bInfo(NULL),linked_bInfo(NULL)
	 {
	    for (int i=0; i<3; ++i)
			{
			_smartSteStart[i]=0;
			_smartSteEnd[i]=0;
			}
	 }
	 
	 
	 //getStencilStart
	 const int * getStencilStart()
	 {
	  return _smartSteStart;
	 }
	 
	 //getSteniclEnd
	 const int * getStencilEnd()
	 {
	  return _smartSteEnd;
	 }
	 
	 //Destructor
	~SmartBoundaryInfo()
	{
		std::cout << "[SmartBoundaries] Calling Destructor (deleting when owing, shallow only when linking)" <<std::endl;

		if (bStatus == updating_owning || bStatus==initialized_owning)
		{
		  assert(owned_bInfo != NULL && linked_bInfo == NULL);
		  delete owned_bInfo;
		  owned_bInfo =NULL;
		  _shallow_clear();
		  
		}
		else if (bStatus == updating_linking || bStatus==initialized_linking)
		{
			assert(owned_bInfo == NULL && linked_bInfo != NULL);
			linked_bInfo =NULL; //just unlink (MRAG will take care of it's boundaries)
			_shallow_clear();
		}
		else
		{
			assert(bStatus == Constructed);
		}
	}
	 
	 
	 //Function Init:
	 template <typename Wavelets, typename Block>
	 void init(Grid<Wavelets,Block>& grid, const int steStart[3], const int steEnd[3])
	 {
	 
	    if (_verb>0) 
		{
			std::cout <<"[smartBoundary]: Initializing .... " <<std::endl;
			std::cout << "requesting stencil: " <<std::endl;
			std::cout << steStart[0] << " " << steStart[1] << " "<< steStart[2] << std::endl;
			std::cout <<  steEnd[0] << " " << steEnd[1] << " "<< steEnd[2] << std::endl;
		}
		

	  
		
		bool reallyOwning=false; //need this due to reinit-checks.
		
		//check OwnerShip.
		BoundaryInfo& tmpInfo=grid.getBoundaryInfo();
		for (int i=0; i<3; ++i)
		{
			if (int(steStart[i])<int(tmpInfo.stencil_start[i]) ||
				int(steEnd[i])>int(tmpInfo.stencil_end[i]))
			{
				reallyOwning=true;
			}			
			
		}
		
		
		//first initialization
		if (bStatus==Constructed)
		{
		
		    assert(linked_bInfo==NULL && owned_bInfo == NULL && bOwnerOfBoundaryInfo==false);
			if(reallyOwning)
			{
    		    if (_verb>0) std::cout << "[SmartBoundaryInfo] First Initialization; Owning " <<std::endl;
				owned_bInfo=grid.createBoundaryInfo(steStart, steEnd);
			}
			else
			{
       		    if (_verb>0) std::cout << "[SmartBoundaryInfo] First Initialization; Linking " <<std::endl;
				linked_bInfo=&grid.getBoundaryInfo();
			}
		}
		
		//reinit or change to different state
		else if(bStatus==initialized_linking || bStatus==updating_linking)
		{
		    
			assert(linked_bInfo!=NULL && owned_bInfo == NULL && bOwnerOfBoundaryInfo==false);
			if (reallyOwning)
			{
				//changed to owning: dont touch the MRAG stuff, refresh stencil, createNewBoundary,updateBase
			    if (_verb>0) std::cout << "[SmartBoundaryInfo] Reinit: Changing from Linking->Owning " <<std::endl;
				linked_bInfo=NULL;
				owned_bInfo=grid.createBoundaryInfo(steStart, steEnd);
			}
			else
			{
			 
				//refresh stencil, dont need to delete, getNewBoundary,updateBase
    		    if (_verb>0) std::cout << "[SmartBoundaryInfo] Reinit: Remains Linking" <<std::endl;
				linked_bInfo=&grid.getBoundaryInfo();
				
			}
		}
		
		//reinit or change to different state
		else if(bStatus==initialized_owning || bStatus==updating_owning)
		{
			assert(linked_bInfo==NULL && owned_bInfo != NULL && bOwnerOfBoundaryInfo==true);
			BoundaryInfo* tmpInfo=owned_bInfo;
			delete tmpInfo;
			if(reallyOwning) //refresh stencil,delete,create new boundary,updateBase
			{
				if (_verb>0) std::cout << "[SmartBoundaryInfo] Reinit: Remains Owning" <<std::endl;
     			owned_bInfo=grid.createBoundaryInfo(steStart,steEnd);
			}
			else  // refresh stencil, delete old pointer, getNewBoundary,updateBase
			{
   				if (_verb>0) std::cout << "[SmartBoundaryInfo] Reinit: Changing from Owning->Linking" <<std::endl;
				owned_bInfo=NULL;
				linked_bInfo=&grid.getBoundaryInfo();
			}
			
		}
		
		
		
	    //Refresh Stencil (in any Case)
		if (reallyOwning)
		{
			for (int i=0; i<3; ++i)
			{
			   stencil_start[i]=char(steStart[i]);
			   _smartSteStart[i]=int(steStart[i]);
			   stencil_end[i]=char(steEnd[i]);
			   _smartSteEnd[i]=int(steEnd[i]);
			}
		}
		else
		{
			 BoundaryInfo& curInfo=grid.getBoundaryInfo();
		     for (int i=0; i<3; ++i)
			 {
				_smartSteStart[i]=int(curInfo.stencil_start[i]);
			   this->stencil_start[i]=curInfo.stencil_start[i];
			 	_smartSteEnd[i]=int(curInfo.stencil_end[i]);
				this->stencil_end[i]=curInfo.stencil_end[i];
			} 
		}
		
		
		//update Ownership
		bOwnerOfBoundaryInfo=reallyOwning;
		//in any case, refreshBase
		if (bOwnerOfBoundaryInfo)
		{ 
			_updateBase(owned_bInfo);
			bStatus=initialized_owning;
		}
		else 
		{
			_updateBase(linked_bInfo);
    		bStatus=initialized_linking;
		}
		
		
		if (_verb>0)
		{						
		 std::cout << "Init finished:" <<std::endl;
         std::cout << "[SmartBoundaries] Debug:, s:" << int(this->stencil_start[0]) << " " <<int(this->stencil_start[1]) << " " << int(this->stencil_start[2]) << std::endl;
 		 std::cout << "[SmartBoundaries] Debug:, e:" << int(this->stencil_end[0]) << " " <<int(this->stencil_end[1]) << " " <<int(this->stencil_end[2]) << std::endl;
		}
	 
	 }
	
	
	 template <typename Wavelets, typename Block>
	 void update(Grid<Wavelets,Block>& grid)
	 { 
		  if (_verb>0) std::cout << "[SmartBoundaries]: starting boundary update...";
		  assert(bStatus!=Constructed); //must have bee initialized 
		  
		  if (bOwnerOfBoundaryInfo) 
		  {
		    assert(bStatus==initialized_owning || bStatus==updating_owning);
			assert(owned_bInfo!=NULL && linked_bInfo==NULL);
			_updateOwning(grid);
		  }
		  else 
		  {	
		    assert(bStatus==initialized_linking || bStatus==updating_linking);
			assert(owned_bInfo==NULL && linked_bInfo!=NULL);
			_updateLinking(grid);
		  }
		 
		 if (_verb>0)
		 {
			 std::cout << "...finished" << std::endl;
			 //std::cout << "[SmartBoundaries] Debug:, s:" << int(this->stencil_start[0]) << " " <<int(this->stencil_start[1]) << " " << int(this->stencil_start[2]) << std::endl;
			 //std::cout << "[SmartBoundaries] Debug:, e:" << int(this->stencil_end[0]) << " " <<int(this->stencil_end[1]) << " " <<int(this->stencil_end[2]) << std::endl;
		 }

	 }
	
	
	}; //End Class SmartBoundaries
	
} //end namespace MRAG