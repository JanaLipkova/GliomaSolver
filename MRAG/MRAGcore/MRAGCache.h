/*
 *  MRAGCache.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 11/4/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <assert.h>

#include <vector>
#include <map>
#undef min
#undef max
#include <queue>
#undef min
#undef max

using namespace std;
namespace MRAG
{

template<typename KeyType, typename CachedType>
class Cache
{
public:
	Cache(void): pcache(NULL), ticket(0), in_cache(),victims(),
		owner(false), valid(false){}

	void setup(vector<CachedType>& cache);
	void setup(const int n);

	CachedType& get(const KeyType key, bool& bCacheHit);

	CachedType& operator[](const KeyType key)
	{
		bool bCacheHit; 
		return get(key, bCacheHit);
	}

	~Cache(void){}

private:

	struct SlotInfo
	{
		bool is_valid;
		unsigned int ticket;
		
		int slot;
		const KeyType* key;
		
		SlotInfo(int ticket_, int slot_, const KeyType * key_=NULL, bool valid=false):
		is_valid(valid),
		ticket(ticket_),
		slot(slot_),
		key(key_)
		{
		}
	
		bool operator < (const SlotInfo& a) const 
		{ 
			if (is_valid) 
				return (!a.is_valid || ticket>a.ticket);
			else
				return false;
		}
	};
	
	void _setup();

	void _rebuildVictims(priority_queue<SlotInfo>& victims) const;

	vector<CachedType> * pcache;
	unsigned int ticket;
	map<KeyType, const int> in_cache;
	priority_queue<SlotInfo> victims;
	bool owner, valid;
};

template< typename KeyType, typename CachedType>
void Cache<KeyType, CachedType>::setup(const int n)
{	
	owner = true;
	
	this->pcache = new vector<CachedType>(n);
	
	_setup();
}

template< typename KeyType, typename CachedType>
void Cache<KeyType, CachedType>::setup(vector<CachedType>& cache)
{
	assert(!valid);
	
	this->pcache = &cache;
	
	_setup();
}

template< typename KeyType, typename CachedType>
void Cache<KeyType, CachedType>::_setup()
{	
	ticket = 0;
	
	in_cache.clear();
	
	victims = priority_queue<SlotInfo>();
	
	const int n = (int)pcache->size();
	for(int i=0; i<n; i++)
		victims.push(SlotInfo(ticket++, i));
	
	valid = true;
}

template< typename KeyType, typename CachedType>
CachedType& Cache<KeyType, CachedType>::get(const KeyType key, bool& bCacheHit)
{
	assert(valid);
	
	CachedType * value = NULL;
	vector<CachedType>& cache = *pcache;
	
	typename map<KeyType, const int>::const_iterator itInCache = in_cache.find(key);
	
	bCacheHit = (itInCache != in_cache.end());
	
	if(bCacheHit)
	{
		//printf("Cache Hit\n");
		assert(itInCache->first == key);
		value = &cache[itInCache->second];
	}
	else
	{
		//printf("Cache Miss\n");
		const SlotInfo& slot_info = victims.top();
		
		value = &cache[slot_info.slot];
		
		if (!slot_info.is_valid)
			*value = CachedType();
		else
			in_cache.erase(*slot_info.key);
		
		in_cache.insert(std::pair<KeyType, const int>(key,slot_info.slot));
		
		victims.push(SlotInfo(ticket++, slot_info.slot, & in_cache.find(key)->first, true));
		victims.pop();
		
		if (ticket == 0)
		{
			_rebuildVictims(victims);
			ticket = (int)victims.size();
		}
	}
	
	return *value;
}

template< typename KeyType, typename CachedType>
void Cache<KeyType, CachedType>::_rebuildVictims(priority_queue<SlotInfo>& v) const
{
	priority_queue<SlotInfo> new_v;
	
	const int n = (int)v.size();
	for(int t=0; t<n; t++)
	{
		SlotInfo slot_info = v.top();
		
		slot_info.ticket = t;
		new_v.push(slot_info);
		
		v.pop();
	}
	
	v = new_v;
}
}