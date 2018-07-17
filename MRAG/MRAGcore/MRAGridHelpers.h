
#include "MRAGrid.h"
#include "MRAGMatrix3D.h"

namespace MRAG
{
	namespace MRAGridHelpers
	{
		inline int _computeDesiredLevel(int nB[3])
		{
			int l=0;
			int ndesired = std::max(nB[0], std::max(nB[1], nB[2]));
			
			for(int n=1; n<ndesired; n*=2) l++;
			
			return l;
		}
		
		template<typename Hierarchy>
		void _createChildren(std::vector<int> & IDs, Hierarchy& h, 
							 GridNode * parent, int child_number, 
							 int n[3], int l, int l_desired)
		{
			int idx[3] = { 0, 0, 0 };
			
			if (parent != NULL)
			{
				idx[0] = parent->index[0]*2 + (child_number&1);
				idx[1] = parent->index[1]*2 + ((child_number>>1)&1);
				idx[2] = parent->index[2]*2 + ((child_number>>2)&1);
			}
			else
				assert(l==0 && child_number==0);
			
			bool bIsOutside =  (idx[0]>= n[0]) || (idx[1]>= n[1]) || (idx[2]>= n[2]);
			if(bIsOutside) return;
			
			if (l<l_desired)
			{
				GridNode * node = new GridNode(true, parent, -1, idx[0], idx[1], idx[2], l);
				
				h[parent].push_back(node);
				h[node].clear();
				
				for(int i=0; i<8; i++)
					_createChildren(IDs, h, node, i, n, l+1, l_desired);
			}
			else
			{
				assert(l==l_desired);
				
				GridNode * node = new GridNode(false, parent, IDs.back(), idx[0], idx[1], idx[2], l_desired);
				
				IDs.pop_back(); 
				
				h[parent].push_back(node);
				h[node].clear();
			}
		}
		
		/*used in compute neighborhood*/
		struct FlattenedBlock
		{
			int index[3];
			int size[3];
			const GridNode * target;
			
			FlattenedBlock(const int reference_start[3], const int reference_end[3],  const GridNode * node):
			target(node)
			{
				index[0] = reference_start[0];	size[0] = reference_end[0] - reference_start[0];
				index[1] = reference_start[1];	size[1] = reference_end[1] - reference_start[1];
				index[2] = reference_start[2];	size[2] = reference_end[2] - reference_start[2];
			}
			
			static bool areAdjacent(const FlattenedBlock& a, const FlattenedBlock& b) 
			{
				const int dX = a.index[0]<b.index[0] ? (b.index[0]-a.index[0]-a.size[0]) : (a.index[0]-b.index[0]-b.size[0]);
				const int dY = a.index[1]<b.index[1] ? (b.index[1]-a.index[1]-a.size[1]) : (a.index[1]-b.index[1]-b.size[1]);
				const int dZ = a.index[2]<b.index[2] ? (b.index[2]-a.index[2]-a.size[2]) : (a.index[2]-b.index[2]-b.size[2]);
				
				return (dX<=0 && dY<=0 && dZ<=0);
			}
			
			static bool find(const GridNode * target , vector<GridNode*>& v)
			{
				const int n = v.size();
				
				for(int i=0; i<n; i++)
					if (v[i] == target) return true;
				
				return false;
			}
			
			template<typename NeighborhoodType>
			static void compareAndAdd(const FlattenedBlock& a, const FlattenedBlock& b,  NeighborhoodType& neighborhood, int &nComparisons, int & nNeighborsConnections)
			{
				GridNode * ta = (GridNode *)a.target;
				GridNode * tb = (GridNode *)b.target;
				nComparisons++;
				if (ta != tb && areAdjacent(a, b) && !find(tb, neighborhood[ta]))
				{
					nNeighborsConnections+=2;
					neighborhood[ta].push_back(tb);
					neighborhood[tb].push_back(ta);
				}
			}
		};
		
	}
}
