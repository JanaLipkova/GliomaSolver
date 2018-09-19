/*
 *  GliomaSolverXcode.cp
 *  GliomaSolverXcode
 *
 *  Created by Lipkova on 19/09/18.
 *  Copyright (c) 2018 Lipkova. All rights reserved.
 *
 */

#include "GliomaSolverXcode.h"
#include "GliomaSolverXcodePriv.h"

CFStringRef GliomaSolverXcodeUUID(void)
{
	CGliomaSolverXcode* theObj = new CGliomaSolverXcode;
	return theObj->UUID();
}

CFStringRef CGliomaSolverXcode::UUID()
{
	return CFSTR("0001020304050607");
}
