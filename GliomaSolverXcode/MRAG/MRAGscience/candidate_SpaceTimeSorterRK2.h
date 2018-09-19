/*
 *  candidate_SpaceTimeSorterRK2.h
 *  MRAG
 *
 *  Created by xxxxxxxxxxxxxxx 4/3/09.
 *  Copyright 2009 Cselab. All rights reserved.
 *
 */

#include "MRAGSpaceTimeSorter.h"

namespace MRAG {
	
	class SpaceTimeSorterRK2 : public SpaceTimeSorter 
		{
		//
	public:
		int estimateNofBlockOperations(int nRecursionsPerLevel, int startLevel) const;
		virtual bool getBlocks(int& level, double& timeStep, double& currentTime, vector<BlockInfo>& blocksInfo,  ETimeInterval& type); 
	};
	
	
	bool SpaceTimeSorterRK2::getBlocks(int& level, double& timeStep, double& currentTime, vector<BlockInfo>& blocksInfo,  ETimeInterval& typeInterval)
	{
		//
		
		//1.
		assert(m_bConnected && m_bInSession);
		assert(m_stackRecursions.empty() == false);
		
		//2.
		const StackItem& item = m_stackRecursions.top();
		const int l = item.level;
		const double t = item.t0;
		const ETimeInterval type = item.type;
		
		//const double dt = m_dt*pow((double)m_nRecursions, -l);
		const double dt = m_dt*pow((double)m_nRecursions, -(l-m_startLevel));
		
		//3.
		level = l;
		timeStep = dt;
		currentTime = t;
		typeInterval = type;
		blocksInfo = (*m_refBlockAtLevel)[l];
		m_nNofPerformedBlockOperations += blocksInfo.size();
		
		if (m_bVerbose)
		{
			switch(type)
			{
				case ETimeInterval_Start:
					printf("L: %d,  t:%f,  dt:%f (%s)\n", l, t, dt, "ETimeInterval_Start");
					break;
				case ETimeInterval_Intermediate1:
					printf("L: %d,  t:%f,  dt:%f (%s)\n", l, t, dt, "ETimeInterval_Intermediate1");
					break;
				case ETimeInterval_End:
					printf("L: %d,  t:%f,  dt:%f (%s)\n", l, t, dt, "ETimeInterval_End");
					break;
				default:
					abort();
			}
		}
		//4.
		m_stackRecursions.pop();
		
		if (type == ETimeInterval_Start)
		{
			m_stackRecursions.push(StackItem(l, t+dt, ETimeInterval_End));
			m_stackRecursions.push(StackItem(l , t+dt, ETimeInterval_Intermediate1));			
			
			if (l < m_maxLevel)
			{
				const double children_dt = dt/m_nRecursions;
				for(int i=m_nRecursions-1; i>=0; i--)
					m_stackRecursions.push(StackItem(l+1 , t+(i+0)*children_dt, ETimeInterval_Start));
			}
		}
		
		//5.
		return (!m_stackRecursions.empty());
	}
	
	int SpaceTimeSorterRK2::estimateNofBlockOperations(int nRecursionsPerLevel, int startLevel) const
	{
		assert(m_bConnected);
		
		stack<StackItem> stackRecursions;
		stackRecursions.push(StackItem(startLevel, 0, ETimeInterval_Start));
		
		const int maxLevel = m_refBlockAtLevel->size() - 1;
		int result = 0;
		
		do
		{
			const StackItem& item = stackRecursions.top();
			const int l = item.level;
			const ETimeInterval type = item.type;
			
			result += (*m_refBlockAtLevel)[l].size();
			
			if (m_bVerbose)
				printf("Estimate: L: %d,  t:FAKE,  dt:FAKE (%s), so far: %d block-ops\n", l, type == ETimeInterval_Start? "Start" : "End", result);
			
			stackRecursions.pop();
			
			if (type == ETimeInterval_Start)
			{
				stackRecursions.push(StackItem(l, 0, ETimeInterval_End));
				stackRecursions.push(StackItem(l , 0, ETimeInterval_Intermediate1));
				
				if (l < maxLevel)
					for(int i=nRecursionsPerLevel-1; i>=0; i--)
						stackRecursions.push(StackItem(l+1 , 0, ETimeInterval_Start));
			}
			
		} while(!stackRecursions.empty());
		
		return result;
	}
}