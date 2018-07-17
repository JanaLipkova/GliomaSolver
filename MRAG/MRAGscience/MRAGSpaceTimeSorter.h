/*
 *  MRAGSpaceTimeSorter.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 5/28/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include <vector>
#include <stack>
using namespace std;

#include "../MRAGcore/MRAGCommon.h"
#include "../MRAGcore/MRAGrid.h"


namespace MRAG
{
    /**
     * Takes care of smaller timesteps in finer grids.
     */
	class SpaceTimeSorter
	{
	public:
		enum ETimeInterval { ETimeInterval_End, ETimeInterval_EndB, ETimeInterval_Start, ETimeInterval_StartB, ETimeInterval_Intermediate1, ETimeInterval_Intermediate2};
		
		template <typename Grid> void connect(Grid& g);
		void disconnect();
		
		
	protected:
		
		struct StackItem
		{
			int level;
			double t0;
			ETimeInterval type;
			StackItem(int l, double t, ETimeInterval type_): level(l), t0(t), type(type_) {}
			StackItem(const StackItem& s):level(s.level), t0(s.t0), type(s.type){}
		};
		
		bool m_bInSession, m_bConnected;
		const bool m_bVerbose;
		int m_nRecursions, m_maxLevel, m_startLevel;
		int m_nNofPerformedBlockOperations, m_nNofExpectedBlockOperations;
		double m_dt, m_t0;
		stack<StackItem> m_stackRecursions;
		
		//not owner
		const vector<vector<BlockInfo> >* m_refBlockAtLevel;
		
		template <typename B, typename W> friend class Grid;
		friend class SpaceTimeSorter_Tester;
		
	private:
		//forbidden
		SpaceTimeSorter(const SpaceTimeSorter&):
		m_bInSession(false), m_bConnected(false),m_bVerbose(false),
		m_refBlockAtLevel(NULL), m_nRecursions(0), m_t0(0), m_dt(0), 
		m_stackRecursions(),m_maxLevel(0){abort();}
		
		SpaceTimeSorter& operator=(const SpaceTimeSorter&){abort(); return *this;}
		
	public:
		
		SpaceTimeSorter():
		m_bInSession(false), m_bConnected(false),m_bVerbose(false),
		m_refBlockAtLevel(NULL), m_nRecursions(0), m_t0(0), m_dt(0), 
		m_stackRecursions(),m_maxLevel(0), m_startLevel(0), m_nNofPerformedBlockOperations(0),m_nNofExpectedBlockOperations(0)
		{
		}
		
		void endSession()
		{
			assert(m_bConnected && m_bInSession);
			assert(m_nNofPerformedBlockOperations == m_nNofExpectedBlockOperations);
			
			m_bInSession = false;
		}
		

		virtual int estimateNofBlockOperations(int nRecursionsPerLevel, int startLevel) const

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
				
				if (false)
					printf("Estimate: L: %d,  t:FAKE,  dt:FAKE (%s), so far: %d block-ops\n", l, type == ETimeInterval_Start? "Start" : "End", result);
				
				stackRecursions.pop();
				
				if (type == ETimeInterval_Start)
				{
					stackRecursions.push(StackItem(l, 0, ETimeInterval_End));
					
					if (l < maxLevel)
						for(int i=nRecursionsPerLevel-1; i>=0; i--)
							stackRecursions.push(StackItem(l+1 , 0, ETimeInterval_Start));
				}
				
			} while(!stackRecursions.empty());
			
			return result;
		}
		
		virtual bool getBlocks(int& level, double& timeStep, double& currentTime, vector<BlockInfo>& blocksInfo, ETimeInterval& typeInterval)
		{
			//1. check status
			//2. retrieve next set of blocks with the correct dt
			//3. fill output
			//4. put next levels to be processed
			//5. return false if not further processing is available
			
			//1.
			assert(m_bConnected && m_bInSession);
			assert(m_stackRecursions.empty() == false);
			
			//2.
			const StackItem& item = m_stackRecursions.top();
			const int l = item.level;
			const double t = item.t0;
			const ETimeInterval type = item.type;
			
			const double dt = m_dt*pow((double)m_nRecursions, -(l-m_startLevel));
			
			//3.
			level = l;
			timeStep = dt;
			currentTime = t;
			typeInterval = type;
			blocksInfo.clear();
			blocksInfo = (*m_refBlockAtLevel)[l];
			m_nNofPerformedBlockOperations += blocksInfo.size();
			
			if (m_bVerbose)
				printf("L: %d,  t:%e,  dt:%e (%s)\n", l, t, dt, type == ETimeInterval_Start? "Start" : "End");
			
			//4.
			m_stackRecursions.pop();
			
			if (type == ETimeInterval_Start)
			{
				m_stackRecursions.push(StackItem(l, t+dt, ETimeInterval_End));
				
				if (l < m_maxLevel)
				{
					const double children_dt = dt/m_nRecursions;
					for(int i=m_nRecursions-1; i>=0; i--)
						m_stackRecursions.push(StackItem(l+1 , t+i*children_dt, ETimeInterval_Start));
				}
			}
			
			//5.
			return (!m_stackRecursions.empty());
		}
		
		void startSession(double timeStep, int nRecursionsPerLevel, double timeStart, int startLevel)
		{
			assert(m_bConnected && !m_bInSession);
			
			m_t0 = timeStart;
			m_dt = timeStep;
			
			m_maxLevel = m_refBlockAtLevel->size() - 1;
			m_startLevel = startLevel;
			
			m_nRecursions = nRecursionsPerLevel;
			m_stackRecursions = stack<StackItem>();
			m_stackRecursions.push(StackItem(startLevel, m_t0, ETimeInterval_Start));
			
			m_nNofPerformedBlockOperations = 0;
			m_nNofExpectedBlockOperations = estimateNofBlockOperations(nRecursionsPerLevel, startLevel);
			
			m_bInSession = true;
		}
	};
	
	template <typename Grid> inline void SpaceTimeSorter::connect(Grid& g)
	{
		assert(!m_bInSession);
		
		m_refBlockAtLevel = &g.m_blockAtLevel;
		
		m_bConnected = true;
	}
	
	inline void SpaceTimeSorter::disconnect()
	{
		assert(!m_bInSession);
		
		m_refBlockAtLevel = NULL;
		
		m_bConnected = false;
	}
}
