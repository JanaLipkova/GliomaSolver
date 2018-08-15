/*
 *  Test.cpp
 *  GliomaXcode
 *
 *  Created by Lipkova on 9/17/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "Test.h"

Test::Test(const int argc, const char ** argv)
{
	var=0;
}

Test::~Test()
{
	std::cout << "----- bye bye ----" << std::endl;
}

void Test::run()
{
    	for(int i = 0; i<1; i++)
		cout << "Required type of simulation is not implemented, please see Glioma_main.cpp for supported types"<< endl;
}