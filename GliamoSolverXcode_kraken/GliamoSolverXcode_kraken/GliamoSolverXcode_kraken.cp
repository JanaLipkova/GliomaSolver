/*
 *  GliamoSolverXcode_kraken.cp
 *  GliamoSolverXcode_kraken
 *
 *  Created by Lipkova on 18/09/18.
 *  Copyright (c) 2018 Lipkova. All rights reserved.
 *
 */

#include "GliamoSolverXcode_kraken.h"
#include "GliamoSolverXcode_krakenPriv.h"

CFStringRef GliamoSolverXcode_krakenUUID(void)
{
	CGliamoSolverXcode_kraken* theObj = new CGliamoSolverXcode_kraken;
	return theObj->UUID();
}

CFStringRef CGliamoSolverXcode_kraken::UUID()
{
	return CFSTR("0001020304050607");
}
