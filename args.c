#include "args.h"
#include <stdio.h>

int argpop(int* pargc, char*** pargv, void* output, arg_type_t type)
{
	//scans the the first argument, and pops it
	//returns 1 if the scan is successful
	fprintf(stderr, "argc == %d\n", *pargc);
	if(*pargc == 0) return 0;
	if(argparse(**pargv, output, type))
	{
		(*pargc)--;
		(*pargv)++;
		return 1;
	}
	return 0;
}

int argparse(char* pargvitem, void* output, arg_type_t type)
{
	char* scanfmt = 0;
	if(type==ARG_INT) scanfmt="%d";
	else if(type==ARG_FLOAT) scanfmt="%f";
	else if(type==ARG_CHAR) scanfmt="%c";
	fprintf(stderr, "What: %s\n", pargvitem);
	int result;
	if(type==ARG_STRING)
	{
		*((char**)output) = pargvitem;
		result = 1;
	}
	else 
	{
		result = sscanf( pargvitem, scanfmt, output);
	}
	return result;
}

int argfind(int* pargc, char***pargv, char* longopt, char* shortopt)
{
	for(int i=0;i<*pargc;i++)
		if(!strcmp((*pargv)[i], longopt) || !strcmp((*pargv)[i], shortopt)) return 1;
	return 0;
}

int argfindval(int* pargc, char***pargv, char* longopt, char* shortopt, arg_type_t type)
{
	for(int i=0;i<*pargc-1;i++)
		if(!strcmp((*pargv)[i], longopt) || !strcmp((*pargv)[i], shortopt)) 
			return argparse(**pargv, (*pargv)[i], type);
	return 0;
}

int gargc;
char** gargv;

int garginit(int argc, char** argv)
{
	gargc = argc;
	gargv = argv;
	return 1;
}

int gargpop(void* output, arg_type_t type)
{
	return argpop(&gargc, &gargv, output, type);
}

int gargfind(char* longopt, char* shortopt)
{
	return argfind(&gargc, &gargv, longopt, shortopt);
}

int gargfindval(char* longopt, char* shortopt, arg_type_t type)
{
	return argfindval(&gargc, &gargv, longopt, shortopt, type);
}

int gargthrow(int n)
{
	gargv+=n;
	gargc-=n;
}
