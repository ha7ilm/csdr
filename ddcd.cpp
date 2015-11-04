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
#include <signal.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <iostream>
#include <vector>
#include <unistd.h>

#define SOFTWARE_NAME "ddcd"
#define MSG_START SOFTWARE_NAME ": "

typedef struct client_s
{
	struct sockaddr_in addr;
	int socket;
	pid_t pid;
	int pipefd[2];
} client_t;

int host_port = 0;
char host_address[100] = "127.0.0.1";
int decimation = 0;

int main(int argc, char* argv[])
{
	int c;
	
	for(;;)
	{
		int option_index = 0;
		static struct option long_options[] = {
		   {"port",       required_argument, 0,  'p' },
		   {"address",    required_argument, 0,  'a' },
		   {"decimation", required_argument, 0,  'd' }
		};
		c = getopt_long(argc, argv, "p:a:d:", long_options, &option_index);
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
		case 0:
		case '?':
		case ':':
		default:
			printf(" 0%o ??\n", c);
		}
	}
	
	if(!decimation) { fprintf(stderr, MSG_START "missing required command line argument, --decimation.\n"); exit(1); }
	if(!host_port) { fprintf(stderr, MSG_START "missing required command line argument, --port.\n"); exit(1); }
	
	struct sockaddr_in addr_host;
    int listen_socket;
	std::vector<client_t*> clients(10);
    listen_socket=socket(AF_INET,SOCK_STREAM,0);
    memset(&addr_host,'0',sizeof(addr_host));
    addr_host.sin_family=AF_INET;
    addr_host.sin_port=htons(host_port);

    if( (addr_host.sin_addr.s_addr=inet_addr(host_address)) == INADDR_NONE) 
	{ fprintf(stderr, MSG_START "invalid host address.\n"); exit(1); }

	if( bind(listen_socket, (struct sockaddr*) &addr_host, sizeof(addr_host)) < 0)
	{ fprintf(stderr, MSG_START "cannot bind() address to the socket.\n"); exit(1); }

	if( listen(listen_socket, 10) == -1)
	{ fprintf(stderr, MSG_START "cannot listen() on socket.\n"); exit(1); }

	for(;;)
	{
		struct sockaddr_in addr_cli;
		socklen_t addr_cli_len;
		int new_socket;

		if( (new_socket = accept(listen_socket, (struct sockaddr*)&addr_cli, &addr_cli_len)) == -1)
		{ 
			fprintf(stderr, MSG_START "cannot accept() a connection.\n"); 
			continue; 
		}

		client_t* new_client = new client_t;
		memcpy(&new_client->addr, &addr_cli, sizeof(new_client->addr));
		new_client->socket = new_socket;
		
		if(new_client->pid = fork())
		{
			//We're the parent
			clients.push_back(new_client);
			printf("client pid: %d\n", new_client->pid);
		}
		else
		{
			//We're the client
			client();
			break;
		}
	}

	return 0;
}

void client()
{
	printf("I'm the client\n");
	for(;;) sleep(1);	
}
