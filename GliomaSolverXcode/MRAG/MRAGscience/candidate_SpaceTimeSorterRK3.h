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
	
	class SpaceTimeSorterRK3 : public SpaceTimeSorter {
		//
	public:
		virtual bool getBlocks(int& level, double& timeStep, double& currentTime, vector<BlockInfo>& blocksInfo,  ETimeInterval& type); 
	};
	
	
	bool SpaceTimeSorterRK3::getBlocks(int& level, double& timeStep, double& currentTime, vector<BlockInfo>& blocksInfo,  ETimeInterval& typeInterval){
		//
		
		//1.
		assert(m_bConnected && m_bInSession);
		assert(m_stackRecursions.empty() == false);
		
		//2.
		const StackItem& item = m_stackRecursions.top();
		const int l = item.level;
		const double t = item.t0;
		const ETimeInterval type = item.type;
		
		const double dt = m_dt*pow((double)m_nRecursions, -l);
		
		//3.
		level = l;
		timeStep = dt;
		currentTime = t;
		typeInterval = type;
		blocksInfo.clear();
		blocksInfo = (*m_refBlockAtLevel)[l];
		
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
				case ETimeInterval_Intermediate2:
					printf("L: %d,  t:%f,  dt:%f (%s)\n", l, t, dt, "ETimeInterval_Intermediate2");
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
			m_stackRecursions.push(StackItem(l , t+0.5*dt, ETimeInterval_Intermediate2));
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
	
    
    
	
	
	
	
	
}