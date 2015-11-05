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


#define SOFTWARE_NAME "ddcd"
#define MSG_START SOFTWARE_NAME ": "

int host_port = 0;
char host_address[100] = "127.0.0.1";
int decimation = 0;
int bufsize = 1024;
int bufsizeall;
char* buf;

int set_nonblocking(int fd)
{
	int flagtmp;
	if((flagtmp = fcntl(fd, F_GETFL))!=-1)
		if((flagtmp = fcntl(fd, F_SETFL, flagtmp|O_NONBLOCK))!=-1)
			return 0;
	return 1;
}

client_t* this_client;

int main(int argc, char* argv[])
{
	int c;
	fd_set select_fds;
	
	for(;;)
	{
		int option_index = 0;
		static struct option long_options[] = {
		   {"port",       required_argument, 0,  'p' },
		   {"address",    required_argument, 0,  'a' },
		   {"decimation", required_argument, 0,  'd' },
		   {"bufsize", required_argument, 	 0,  'b' }
		};
		c = getopt_long(argc, argv, "p:a:d:b:", long_options, &option_index);
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
		case 0:
		case '?':
		case ':':
		default:
			printf(" 0%o ??\n", c);
		}
	}
	
	if(!decimation) error_exit(MSG_START "missing required command line argument, --decimation.\n");
	if(!host_port) error_exit(MSG_START "missing required command line argument, --port.\n");
	
	struct sockaddr_in addr_host;
    int listen_socket;
	std::vector<client_t*> clients;
	clients.reserve(100);
    listen_socket=socket(AF_INET,SOCK_STREAM,0);

	int sockopt = 1;
	if( setsockopt(listen_socket, SOL_SOCKET, SO_REUSEADDR, (char *)&sockopt, sizeof(sockopt)) == -1 )
		error_exit(MSG_START "cannot set SO_REUSEADDR.\n");  //the best description on SO_REUSEADDR ever: http://stackoverflow.com/a/14388707/3182453

    memset(&addr_host,'0',sizeof(addr_host));
    addr_host.sin_family=AF_INET;
    addr_host.sin_port=htons(host_port);
	addr_host.sin_addr.s_addr = INADDR_ANY;

    if( (addr_host.sin_addr.s_addr=inet_addr(host_address)) == INADDR_NONE ) 
		error_exit(MSG_START "invalid host address.\n");

	if( bind(listen_socket, (struct sockaddr*) &addr_host, sizeof(addr_host)) < 0 )
		error_exit(MSG_START "cannot bind() address to the socket.\n");

	if( listen(listen_socket, 10) == -1 )
		error_exit(MSG_START "cannot listen() on socket.\n");

	struct sockaddr_in addr_cli;
	socklen_t addr_cli_len = sizeof(addr_cli);
	int new_socket;

	//The server will wait on these sockets later...


	//Set stdin and listen_socket to non-blocking
	if(set_nonblocking(STDIN_FILENO) || set_nonblocking(listen_socket))
		error_exit(MSG_START "cannot set_nonblocking().\n");

	bufsizeall = bufsize*sizeof(char);
	buf = (char*)malloc(bufsizeall);

	FD_ZERO(&select_fds);
	FD_SET(listen_socket, &select_fds);
	FD_SET(STDIN_FILENO, &select_fds);
	int highfd = ((listen_socket>STDIN_FILENO)?listen_socket:STDIN_FILENO)  + 1;

	for(;;)
	{
		//Let's wait until there is any new data to read, or any new connection!
		select(highfd, &select_fds, NULL, NULL, NULL);

		//Is there a new client connection?
		if( (new_socket = accept(listen_socket, (struct sockaddr*)&addr_cli, &addr_cli_len)) != -1)
		{ 
			this_client = new client_t;
			memcpy(&this_client->addr, &addr_cli, sizeof(this_client->addr));
			this_client->socket = new_socket;
			if(pipe(this_client->pipefd) == -1)
			{ 
				perror(MSG_START "cannot open new pipe() for the client.\n");
				continue;
			}
			if(this_client->pid = fork())
			{
				//We're the parent
				set_nonblocking(this_client->pipefd[1]);
				clients.push_back(this_client);
				printf("client pid: %d\n", this_client->pid);
			}
			else
			{
				//We're the client
				client();
				return 1;
			}
		}

		int retval = read(STDIN_FILENO, buf, bufsizeall);
		if(retval==0)
		{
			//end of input stream, close clients and exit
		}
		else if(retval != -1)
		{
			for (int i=0;i<clients.size(); i++)
			{
				if(write(clients[i]->pipefd[1], buf, retval)==-1)
					print_client(clients[i], "lost buffer, failed to write pipe");
			}
		}
		//TODO: at the end, server closes pipefd[1] for client
		//close(this_client->pipefd[1]);
	}

	return 0;
}

void print_client(client_t* client, const char* what)
{
	fprintf(stderr,MSG_START " (client %s:%d) %s\n", inet_ntoa(client->addr.sin_addr), client->addr.sin_port, what);
}

void client_cleanup()
{
	close(this_client->pipefd[0]);
}

void client()
{
	printf("I'm the client\n");
	for(;;)
	{
		read(this_client->pipefd[0],buf,bufsizeall);
		send(this_client->socket,buf,bufsizeall,0);
	}	
}

void error_exit(const char* why)
{
	perror(why);
	exit(1);
}
