/*
This software is part of libcsdr, a set of simple DSP routines for 
Software Defined Radio.

Copyright (c) 2014, Andras Retzler <randras@sdr.hu>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the copyright holder nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ANDRAS RETZLER BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "ddcd.h"

int host_port = 0;
char host_address[100] = "127.0.0.1";
int decimation = 0;
float transition_bw = 0.05;
int bufsize = 1024;
char ddc_method_str[100] = "td";
ddc_method_t ddc_method;

int main(int argc, char* argv[])
{
	int c;
	for(;;)
	{
		int option_index = 0;
		static struct option long_options[] = {
		   {"port",       required_argument, 0,  'p' },
		   {"address",    required_argument, 0,  'a' },
		   {"decimation", required_argument, 0,  'd' },
		   {"bufsize", 	  required_argument, 0,  'b' },
	       {"method", 	  required_argument, 0,  'm' },
	       {"transition", required_argument, 0,  't' }
		};
		c = getopt_long(argc, argv, "p:a:d:b:m:t:", long_options, &option_index);
		if(c==-1) break;
		switch (c) 
		{
		case 'a':
			host_address[100-1]=0;
			strncpy(host_address,optarg,100-1);
			break;
		case 'p':
			host_port=atoi(optarg);
			break;
		case 'd':
			decimation=atoi(optarg);
			break;
		case 'b':
			bufsize=atoi(optarg);
			break;
		case 'm':
			ddc_method_str[100-1]=0;
			strncpy(ddc_method_str,optarg,100-1);
			break;
		case 't':
			sscanf(optarg,"%g",&transition_bw);
			break;
		case 0:
		case '?':
		case ':':
		default:;
			print_exit(MSG_START "error in getopt_long()\n");
		}
	}
	
	if(!decimation) print_exit(MSG_START "missing required command line argument, --decimation.\n");
	if(!host_port) print_exit(MSG_START "missing required command line argument, --port.\n");
	if(decimation<0) print_exit(MSG_START "invalid value for --decimation (should be >0).\n");
	if(decimation==1) fprintf(stderr, MSG_START "decimation = 1, just copying raw samples.\n");
	if(transition_bw<0||transition_bw>0.5) print_exit(MSG_START "invalid value for --transition (should be between 0 and 0.5).\n");
	
	if(decimation==1); //don't do anything then
	else if(!strcmp(ddc_method_str,"td")) 
	{
		ddc_method = M_TD; 
		fprintf(stderr, MSG_START "method is M_TD (default).\n");
	}
	else if (!strcmp(ddc_method_str,"fastddc")) 
	{
		ddc_method = M_FASTDDC; 
		fprintf(stderr, MSG_START "method is M_FASTDDC.\n");
	}
	else print_exit(MSG_START "invalid parameter given to --method.\n");


}

void print_exit(const char* why)
{
	fprintf(stderr, "%s", why);
	exit(1);
}
