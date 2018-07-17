/*
 *  MRAGNode.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 4/24/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once

namespace MRAG
{
	struct GridNode
	{
		bool isEmpty;
		bool shouldBeCompressed;
		bool shouldBeRefined;
		GridNode * parent;
		int blockID;
		int index[3];
		int level;

		GridNode(bool isEmpty_, GridNode* parent_, int blockID_, int idX, int idY, int idZ, int level_):
			isEmpty(isEmpty_),
			parent(parent_),
			blockID(blockID_),
			level(level_),
			shouldBeCompressed(false),
			shouldBeRefined(false)
		{
			index[0] = idX;
			index[1] = idY;
			index[2] = idZ;
		}
	};
}

