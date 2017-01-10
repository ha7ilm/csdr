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

#include "nmux.h"

int host_port = 0;
char host_address[100] = "127.0.0.1";
int thread_cntr = 0;

//CLI parameters
int bufsize = 1024; //! currently unused
int bufcnt = 1024;

char** global_argv;
int global_argc;
tsmpool* pool;

void sig_handler(int signo)
{
	fprintf(stderr, MSG_START "signal %d caught, exiting...\n", signo);
	fflush(stderr);
	exit(0);
}

int main(int argc, char* argv[])
{
	global_argv = argv;
	global_argc = argc;
	int c;
	for(;;)
	{
		int option_index = 0;
		static struct option long_options[] = {
		   {"port",       required_argument, 0,  'p' },
		   {"address",    required_argument, 0,  'a' },
		   {"bufsize", 	  required_argument, 0,  'b' },
		   {"bufcnt", 	  required_argument, 0,  'n' }
		};
		c = getopt_long(argc, argv, "p:a:b:n:", long_options, &option_index);
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
		case 'b':
			bufsize=atoi(optarg);
			break;
		case 'n':
			bufcnt=atoi(optarg);
			break;
		case 0:
		case '?':
		case ':':
		default:;
			print_exit(MSG_START "error in getopt_long()\n");
		}
	}

	if(!host_port) print_exit(MSG_START "missing required command line argument, --port.\n");
	if(bufsize<0) print_exit(MSG_START "invalid value for --bufsize (should be >0)\n");
	if(bufcnt<0) print_exit(MSG_START "invalid value for --bufcnt (should be >0)\n");

	//set signals
	struct sigaction sa;
	memset(&sa, 0, sizeof(sa));
	sa.sa_handler = sig_handler;
	sigaction(SIGTERM, &sa, NULL);
	sigaction(SIGKILL, &sa, NULL);
	sigaction(SIGQUIT, &sa, NULL);
	sigaction(SIGINT, &sa, NULL);
	sigaction(SIGHUP, &sa, NULL);

	struct sockaddr_in addr_host;
    int listen_socket;
	std::vector<client_t*> clients;
	clients.reserve(100);
    listen_socket=socket(AF_INET,SOCK_STREAM,0);

	int sockopt = 1;
	if( setsockopt(listen_socket, SOL_SOCKET, SO_REUSEADDR, (char *)&sockopt, sizeof(sockopt)) == -1 )
		error_exit(MSG_START "cannot set SO_REUSEADDR");  //the best description on SO_REUSEADDR ever: http://stackoverflow.com/a/14388707/3182453

	memset(&addr_host,'0',sizeof(addr_host));
    addr_host.sin_family = AF_INET;
    addr_host.sin_port = htons(host_port);
	addr_host.sin_addr.s_addr = INADDR_ANY;

    if( (addr_host.sin_addr.s_addr=inet_addr(host_address)) == INADDR_NONE )
		error_exit(MSG_START "invalid host address");

	if( bind(listen_socket, (struct sockaddr*) &addr_host, sizeof(addr_host)) < 0 )
		error_exit(MSG_START "cannot bind() address to the socket");

	if( listen(listen_socket, 10) == -1 )
		error_exit(MSG_START "cannot listen() on socket");

	fprintf(stderr,MSG_START "listening on %s:%d\n", inet_ntoa(addr_host.sin_addr), host_port);

	struct sockaddr_in addr_cli;
	socklen_t addr_cli_len = sizeof(addr_cli);
	int new_socket;

	int highfd = 0;
	FD_ZERO(&select_fds);
	FD_SET(listen_socket, &select_fds);
	maxfd(&highfd, listen_socket);
	FD_SET(input_fd, &select_fds);
	maxfd(&highfd, input_fd);

	//Set stdin and listen_socket to non-blocking
	if(set_nonblocking(input_fd) || set_nonblocking(listen_socket))
		error_exit(MSG_START "cannot set_nonblocking()");

	//Create tsmpool
	pool = new tsmpool(bufsize, bufcnt);
	if(!pool->ok) print_exit(MSG_START "tsmpool failed to initialize\n");

	unsigned char* current_write_buffer = pool->get_write_buffer();
	int index_in_current_write_buffer = 0;

	pthread_cond_t* wait_condition = new pthread_cond_t;
	if(!pthread_cond_init(wait_condition, NULL)) 
		print_exit(MSG_START "pthread_cond_init failed"); //cond_attrs is ignored by Linux

	for(;;)
	{
		//Let's wait until there is any new data to read, or any new connection!
		select(highfd, &select_fds, NULL, NULL, NULL);

		//Is there a new client connection?
		if( (new_socket = accept(listen_socket, (struct sockaddr*)&addr_cli, &addr_cli_len)) != -1)
		{
			clients_close_all_finished();

			//We're the parent, let's create a new client and initialize it
			client_t* new_client = new client_t;
			new_client->error = 0;
			memcpy(&new_client->addr, &addr_cli, sizeof(new_client->addr));
			new_client->socket = new_socket;
			new_client->status = CS_CREATED;
			new_client->tsmthread = pool->register_thread();
			new_client->lpool = pool;
			pthread_mutex_init(&new_client->wait_mutex, NULL);
			new_client->wait_condition = wait_condition;
			new_client->sleeping = 0;
			if(pthread_create(&new_client->thread, NULL, client_thread , (void*)&new_client)<0)
			{
				clients.push_back(new_client);
				fprintf(stderr, MSG_START "pthread_create() done, clients now: %d\n", clients.size());
			}
			else
			{
				fprintf(stderr, MSG_START "pthread_create() failed.\n");
			}
		}

		if(index_in_current_write_buffer >= bufsize)
		{
			current_write_buffer = pool->get_write_buffer();
			pthread_cond_broadcast(wait_condition);
			index_in_current_write_buffer = 0;
		}
		int retval = read(input_fd, current_write_buffer + index_in_current_write_buffer, bufsize - index_in_current_write_buffer);
		if(retval>0)
		{
			index_in_current_write_buffer += retval;
		}
		else if(retval==0)
		{
			//!end of input stream, close clients and exit
			print_exit(MSG_START "end of input, exiting.\n")
		}
	}
}

#if 0
for (int i=0; i<clients.size(); i++)
{
	if(write(clients[i]->pipefd[1], buf, retval)==-1)
	{

		if(!clients[i]->error)
		{
			print_client(clients[i], "lost buffer, failed to write pipe.");
			clients[i]->error=1;
		}
		//fprintf(stderr, MSG_START "errno is %d\n", errno); //usually 11
		//int wpstatus;
		//int wpresult = waitpid(clients[i]->pid, &wpstatus, WNOHANG);
		//fprintf(stderr, MSG_START "pid is %d\n",clients[i]->pid);
		//perror("somethings wrong");
		//if(wpresult == -1) print_client(clients[i], "error while waitpid()!");
		//else if(wpresult == 0)
		waitpid(clients[i]->pid, NULL, WNOHANG);
		if(!proc_exists(clients[i]->pid))
		{
			//Client exited!
			print_client(clients[i], "closing client from main process.");
			close(clients[i]->pipefd[1]);
			close(clients[i]->socket);
			delete clients[i];
			clients.erase(clients.begin()+i);
			fprintf(stderr, MSG_START "done closing client from main process.\n");
		}
	}
	else  { if(clients[i]->error) print_client(clients[i], "pipe okay again."); clients[i]->error=0; }
}
}
//TODO: at the end, server closes pipefd[1] for client
#endif

void client_erase(client_t* client)
{
	pthread_mutex_destroy(client->wait_mutex);
	pthread_cond_destroy(client->wait_condition);
	pool->
}

void clients_close_all_finished()
{
	for(int i=0;i<clients.size();i++)
	{
		if(clients[i]->status == CS_THREAD_FINISHED)
		{
			client_erase(clients[i]);
			clients.erase(i);
		}
	}
}

void* client_thread (void* param)
{
	client_t* this_client = (client_t*)param;
	this_client->status = CS_THREAD_RUNNING;
	int retval;
	tsmpool* lpool = this_client->lpool;
	for(;;)
	{
		//wait until there is any data in the tsmpool for me (wait for the server process to wake me up)
		for(;;)
		{
			char* pool_read_buffer = (char*)lpool->get_read_buffer();
			if(pool_read_buffer) break;
			pthread_mutex_lock(&this_client->wait_mutex);
			this_client->sleeping = 1;
			pthread_cond_wait(this_client->wait_condition, &this_client->wait_mutex);
		}

		//wait for the socket to be available for write
		//read data from server global tsmpool
		//write data to client socket
		//have an exit condition


	}
	this_client->status = CS_THREAD_FINISHED;
	pthread_exit(NULL);
	return NULL;
}

void error_exit(const char* why)
{
	perror(why); //do we need a \n at the end of (why)?
	exit(1);
}

void print_exit(const char* why)
{
	fprintf(stderr, "%s", why);
	exit(1);
}

void maxfd(int* maxfd, int fd)
{
	if(fd>=*maxfd) *maxfd=fd+1;
}
