/*
 *  MRAG_STDTestL1_BlockCollection.h
 *  MRAG
 *
 *  Created by Thycall Ya'acov on 7/21/08.
 *  Copyright 2008 MichaelBergdorf. All rights reserved.
 *
 */

#ifndef MRAG_STDTestL1_BlockCollection_
#define MRAG_STDTestL1_BlockCollection_
#include "MRAG_STDTestL1.h"
#include "MRAGBlockCollection.h"
#include <iostream>

namespace MRAG {
	
	
	template <typename Wavelets, typename Block>
	class MRAG_STDTestL1_BlockCollection : public MRAG_STDTestL1<Wavelets,Block> {
	protected:
		void create(BlockCollection<Block> &bc, set<int> &blockIDSet, unsigned int howmany)
		{
			vector<int> newBlocks = bc.create(howmany);
			vector<int>::iterator it;
			for(vector<int>::iterator it = newBlocks.begin();it!=newBlocks.end();++it){
				blockIDSet.insert(*it);
			}
		//	std::cout << "inserted " << newBlocks.size() << " blocks" << std::endl;
		}
		
		void remove(BlockCollection<Block> &bc, std::set<int> &blockIDSet, unsigned int howmany)
		{
			howmany = std::min((int)howmany, (int)blockIDSet.size());
			
			for(int iRemove=0; iRemove<howmany; iRemove++)
			{
				std::set<int>::const_iterator it = blockIDSet.begin();
				
				const int index = rand() % blockIDSet.size();
				for(int i=0; i<index; i++) it++;
				
				assert(it!=blockIDSet.end());
				
				bc.erase(*it);
				blockIDSet.erase(*it);
			}
			
		//	std::cout << "removed " << realHowmany << " blocks" << std::endl;
		}
			
			
		bool run()
		{
			try	{
				/* 
				 
				 STANDARD OUTPUT

				 MRAG::MRAG_STDTestL1_BlockCollection
				 1 blocks, memorySize = 0.000240326 vec.size = 1 while nChunkSize is 3
				 2 blocks, memorySize = 0 vec.size = 2 while nChunkSize is 3
				 4 blocks, memorySize = 0.000480652 vec.size = 4 while nChunkSize is 3
				 8 blocks, memorySize = 0.000480652 vec.size = 8 while nChunkSize is 3
				 16 blocks, memorySize = 0.00144196 vec.size = 16 while nChunkSize is 3
				 32 blocks, memorySize = 0.00240326 vec.size = 32 while nChunkSize is 3
				 64 blocks, memorySize = 0.00528717 vec.size = 64 while nChunkSize is 3
				 128 blocks, memorySize = 0.0100937 vec.size = 128 while nChunkSize is 3
				 256 blocks, memorySize = 0.020668 vec.size = 256 while nChunkSize is 3
				 512 blocks, memorySize = 0.0408554 vec.size = 512 while nChunkSize is 3
				 1024 blocks, memorySize = 0.0821915 vec.size = 1024 while nChunkSize is 3
				 2048 blocks, memorySize = 0.163902 vec.size = 2048 while nChunkSize is 3
				 4096 blocks, memorySize = 0.328285 vec.size = 4096 while nChunkSize is 3
				 8192 blocks, memorySize = 0.65609 vec.size = 8192 while nChunkSize is 3
				 16384 blocks, memorySize = 1.31266 vec.size = 16384 while nChunkSize is 3
				 32768 blocks, memorySize = 2.62484 vec.size = 32768 while nChunkSize is 3
				 65536 blocks, memorySize = 5.25016 vec.size = 65536 while nChunkSize is 3
				 131072 blocks, memorySize = 10.4998 vec.size = 131072 while nChunkSize is 3
				 262144 blocks, memorySize = 21.0002 vec.size = 262144 while nChunkSize is 3
				 524288 blocks, memorySize = 41.9998 vec.size = 524288 while nChunkSize is 3
				 
				 */
				
				BlockCollection<Block> blockCollection;
				const unsigned int number = 10;
				unsigned int nBlocks = 1 << number;
				set<int> blockSet;
				vector<int>::iterator it;
				vector<int> blocks = blockCollection.create(nBlocks);
				for(it=blocks.begin();it!=blocks.end();++it) blockSet.insert(*it);

				/* runn our tests now */
				for(int i=0;i<200;++i){
					create(blockCollection,blockSet,100);
					//int i;
					printf("Collection size %f MB\n", blockCollection.getMemorySize());
					//cin>>i;
					remove(blockCollection,blockSet,100);
					
					
					//cin>>i;
				}
				std::cout << "Test complete" <<std::endl;
				
				
				blockCollection.clear();
				
				blocks.clear();
				
				//exit(0);
				blockSet.clear();
				
				
			}
			catch(...)
			{
				return false;
			}
			
			return true;
		}
	public:
		static void runTests(){
			MRAG_STDTestL1_BlockCollection<Wavelets, Block> test;
			assert(test.run());
		}
		
	}; 
	
}


#endif