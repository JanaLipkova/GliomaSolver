/*
 *  Test.h
 *  GliomaXcode
 *
 *  Created by Lipkova on 9/17/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#pragma once
#include "Glioma.h"

using namespace MRAG;

class Test: public Glioma
{
	int var;
public:
	Test(const int argc, const char ** argv);
	~Test();
	void run();
};
	
