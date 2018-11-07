/*
 *  ArgumentParser.h
 *
 *	Simplified version of the argument parser
 *
 *  Created by Christian Conti on 9/26/2012
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */


//========================================================================================
//
//	Argument Parser
//		this two classes (Value, ArgumentParser) are used to parse parameters
//			from command line
//		usage from command line:
//			after calling the program list the argument as shown in the example below:
//				"./program -parameter1 value1 -parameter2 value2"
//		usage from within the program:
//			(i)  ArgumentParser parser(argc, argv)
//			(ii) parser("-ParameterName").asType(defaultValue)
//				 where Type is Bool, Double, Int or String
//
//========================================================================================

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <vector>
#include <string>
#include <cassert>
#include <sys/time.h>

using namespace std;


class Value
{
private:
	string content;
	
public:
	
	Value() : content("") {}
	
	Value(string content_) : content(content_) {}
	
	double asDouble(double def=0) const
	{
		if (content == "") return def;
		return (double) atof(content.c_str());
	}
	
	int asInt(int def=0) const
	{
		if (content == "") return def;
		return atoi(content.c_str());
	}
	
	bool asBool(bool def=false) const 
	{ 
		if (content == "") return def; 
		if (content == "0") return false; 
		if (content == "false") return false; 
		
		return true;
	}
	
	string asString(string def="") const 
	{ 
		if (content == "") return def;
		
		return content; 
	}
};

class ArgumentParser
{
private:
	
	map<string,Value> mapArguments;
	
	const int iArgC;
	const char** vArgV;
	
public:
	
	Value operator()(const string arg)
	{
		printf("%s is %s\n", arg.data(), mapArguments[arg].asString().data());
		return mapArguments[arg];
	}
	
	ArgumentParser(const int argc, const char ** argv) : mapArguments(), iArgC(argc), vArgV(argv)
	{
		for (int i=1; i<argc; i++)
			if (argv[i][0] == '-')
			{
				string values = "";
				int itemCount = 0;
				
				for (int j=i+1; j<argc; j++)
					if (argv[j][0] == '-')
						break;
					else
					{
						if (strcmp(values.c_str(), ""))
							values += ' ';
						
						values += argv[j];
						itemCount++;
					}
				
				if (itemCount == 0)
					values += '1';
				mapArguments[argv[i]] = Value(values);
				i += itemCount;
			}
	}
};