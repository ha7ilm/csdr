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

char help_text[]="nmux is a TCP stream multiplexer. It reads data from the standard input, and sends it to each client connected through TCP sockets. Available command line options are:\n"
"\t--port (-p), --address (-a): TCP port and address to listen.\n"
"\t--bufsize (-b), --bufcnt (-n): Internal buffer size and count.\n"
"\t--help (-h): Show this message.\n";

int host_port = 0;
char host_address[100] = "127.0.0.1";
int thread_cntr = 0;

//CLI parameters
int bufsize = 1024; 
int bufcnt = 1024;

char** global_argv;
int global_argc;
tsmpool* pool;

pthread_cond_t wait_condition;
pthread_mutex_t wait_mutex;

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
	int no_options = 1;
	for(;;)
	{
		int option_index = 0;
		static struct option long_options[] = {
		   {"port",       required_argument, 0,  'p' },
		   {"address",    required_argument, 0,  'a' },
		   {"bufsize", 	  required_argument, 0,  'b' },
		   {"bufcnt", 	  required_argument, 0,  'n' },
		   {"help", 	  no_argument, 		 0,  'h' },
		   {0,			  0,                 0,  0   }
		};
		c = getopt_long(argc, argv, "p:a:b:n:h", long_options, &option_index);
		if(c==-1) break;
		no_options = 0;
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
		case 'h':
			print_exit(help_text);
			break;
		case 0:
		case '?':
		case ':':
		default:
			print_exit(MSG_START "error in getopt_long()\n");
		}
	}

	if(no_options) print_exit(help_text);
	if(!host_port) print_exit(MSG_START "missing required command line argument, --port.\n");
	if(bufsize<=0) print_exit(MSG_START "invalid value for --bufsize (should be >0)\n");
	if(bufcnt<=0) print_exit(MSG_START "invalid value for --bufcnt (should be >0)\n");

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

	fprintf(stderr, MSG_START "listening on %s:%d\n", inet_ntoa(addr_host.sin_addr), host_port);

	struct sockaddr_in addr_cli;
	socklen_t addr_cli_len = sizeof(addr_cli);
	int new_socket;

	int highfd = 0;
	maxfd(&highfd, listen_socket);
	maxfd(&highfd, STDIN_FILENO);
	
	fd_set select_fds;

	//Set stdin and listen_socket to non-blocking
	if(set_nonblocking(STDIN_FILENO) || set_nonblocking(listen_socket))
		error_exit(MSG_START "cannot set_nonblocking()");

	//Create tsmpool
	pool = new tsmpool(bufsize, bufcnt);
	if(!pool->is_ok()) print_exit(MSG_START "tsmpool failed to initialize\n");

	unsigned char* current_write_buffer = (unsigned char*)pool->get_write_buffer();
	int index_in_current_write_buffer = 0;

	//Create wait condition: client threads waiting for input data from the main thread will be
	//	waiting on this condition. They will be woken up with pthread_cond_broadcast() if new
	//	data arrives.
	if(pthread_cond_init(&wait_condition, NULL)) 
		print_exit(MSG_START "pthread_cond_init failed"); //cond_attrs is ignored by Linux
	
	if(pthread_mutex_init(&wait_mutex, NULL))
		print_exit(MSG_START "pthread_mutex_t failed"); //cond_attrs is ignored by Linux

	for(;;)
	{
		FD_ZERO(&select_fds);
		FD_SET(listen_socket, &select_fds);
		FD_SET(STDIN_FILENO, &select_fds);

		if(NMUX_DEBUG) fprintf(stderr, "mainfor: selecting...");
		//Let's wait until there is any new data to read, or any new connection!
		int select_ret = select(highfd, &select_fds, NULL, NULL, NULL);
		if(NMUX_DEBUG) fprintf(stderr, "selected.\n");
		if(select_ret == -1) error_exit("mainfor select() error");

		//Is there a new client connection?
		if( FD_ISSET(listen_socket, &select_fds) && ((new_socket = accept(listen_socket, (struct sockaddr*)&addr_cli, &addr_cli_len)) != -1) )
		{
			if(NMUX_DEBUG) 
			{
				fprintf(stderr, "\x1b[1m\x1b[33mmainfor: clients before closing: ");
				for(int i=0;i<clients.size();i++) fprintf(stderr, "%p ", clients[i]);
				fprintf(stderr, "\x1b[0m\n");
			}
			if(NMUX_DEBUG) fprintf(stderr, "mainfor: accepted (socket = %d).\n", new_socket);
			//Close all finished clients
			for(int i=0;i<clients.size();i++)
			{
				if(clients[i]->status == CS_THREAD_FINISHED)
				{
                    if(pthread_detach(clients[i]->thread)!=0)
                    {
                        fprintf(stderr,"nmux pthread_detach failed for client %d\n", i);
                        continue;
                    }

					if(NMUX_DEBUG) fprintf(stderr, "mainfor: client removed: %d\n", i);
					//client destructor
					pool->remove_thread(clients[i]->tsmthread);
					clients.erase(clients.begin()+i);
					i--;
				}
			}
			if(NMUX_DEBUG) 
			{
				fprintf(stderr, "\x1b[1m\x1b[33mmainfor: clients after closing: ");
				for(int i=0;i<clients.size();i++) fprintf(stderr, "%p ", clients[i]);
				fprintf(stderr, "\x1b[0m\n");
			}

			//We're the parent, let's create a new client and initialize it
			client_t* new_client = new client_t;
			new_client->error = 0;
			memcpy(&new_client->addr, &addr_cli, sizeof(struct sockaddr_in));
			new_client->socket = new_socket;
			new_client->status = CS_CREATED;
			new_client->tsmthread = pool->register_thread();
			new_client->lpool = pool;
			new_client->sleeping = 0;
			if(pthread_create(&new_client->thread, NULL, client_thread, (void*)new_client)==0)
			{
				clients.push_back(new_client);
				fprintf(stderr, MSG_START "pthread_create() done, clients now: %d\n", (int)clients.size());
			}
			else
			{
				fprintf(stderr, MSG_START "pthread_create() failed.\n");
				pool->remove_thread(new_client->tsmthread);
				delete new_client;
			}
		}

		if( FD_ISSET(STDIN_FILENO, &select_fds) )
		{
			if(index_in_current_write_buffer >= bufsize)
			{
				if(NMUX_DEBUG) fprintf(stderr, "mainfor: gwbing...");
				current_write_buffer = (unsigned char*)pool->get_write_buffer();
				if(NMUX_DEBUG) fprintf(stderr, "gwbed.\nmainfor: cond broadcasting...");
				pthread_mutex_lock(&wait_mutex);
				pthread_cond_broadcast(&wait_condition); 
				pthread_mutex_unlock(&wait_mutex);
				if(NMUX_DEBUG) fprintf(stderr, "cond broadcasted.\n");
					//Shouldn't we do it after we put data in?
					//	No, on get_write_buffer() actually the previous buffer is getting available 
					//	for read for threads that wait for new data (wait on global pthead mutex 
					//	wait_condition). 
				index_in_current_write_buffer = 0;
			}

			if(NMUX_DEBUG) fprintf(stderr, "mainfor: reading...\n");
			int read_ret = read(STDIN_FILENO, current_write_buffer + index_in_current_write_buffer, bufsize - index_in_current_write_buffer);
			if(NMUX_DEBUG) fprintf(stderr, "read %d\n", read_ret);
			if(read_ret>0)
			{
				index_in_current_write_buffer += read_ret;
			}
			else if(read_ret==0)
			{
				//End of input stream, close clients and exit
				print_exit(MSG_START "(main thread/for) end input stream, exiting.\n");
			}
			else if(read_ret==-1)
			{
				if(errno == EAGAIN) { if(NMUX_DEBUG) fprintf(stderr, "mainfor: read EAGAIN\n"); /* seems like select would block forever, so we just read again */ }
				else error_exit(MSG_START "(main thread/for) error in read(), exiting.\n");
			}
		}
	}
}

void* client_thread (void* param)
{
	fprintf(stderr, "client %p: started!\n", param);
	client_t* this_client = (client_t*)param;
	this_client->status = CS_THREAD_RUNNING;
	int retval;
	tsmpool* lpool = this_client->lpool;
	if(NMUX_DEBUG) fprintf(stderr, "client %p: socket = %d!\n", param, this_client->socket);

	if(NMUX_DEBUG) fprintf(stderr, "client %p: poll init...", param);
	struct pollfd pollfds[1];
	pollfds[0].fd = this_client->socket;
	pollfds[0].events = POLLOUT;
	pollfds[0].revents = 0;
	if(NMUX_DEBUG) fprintf(stderr, "client poll inited.\n");

	//Set this_client->socket to non-blocking
	if(set_nonblocking(this_client->socket))
		error_exit(MSG_START "cannot set_nonblocking() on this_client->socket");

	int client_buffer_index = 0;
	int client_goto_source = 0;
	char* pool_read_buffer = NULL;

	for(;;)
	{
		//Wait until there is any data to send.
		//  If I haven't sent all the data from my last buffer, don't wait.
		//	(Wait for the server process to wake me up.)
		while(!pool_read_buffer || client_buffer_index >= lpool->size)
		{
			if(NMUX_DEBUG) fprintf(stderr, "client %p: trying to grb\n", param);
			pool_read_buffer = (char*)lpool->get_read_buffer(this_client->tsmthread);
			if(pool_read_buffer) { client_buffer_index = 0; break; }
			if(NMUX_DEBUG) fprintf(stderr, "client %p: cond_waiting for more data\n", param);
			pthread_mutex_lock(&wait_mutex);
			this_client->sleeping = 1;
			pthread_cond_wait(&wait_condition, &wait_mutex);
			pthread_mutex_unlock(&wait_mutex);
		}

		//Wait for the socket to be available for write.
		if(NMUX_DEBUG) fprintf(stderr, "client %p: polling for socket write...", param);
		int ret = poll(pollfds, 1, -1);
		if(NMUX_DEBUG) fprintf(stderr, "client polled for socket write.\n");
		if(ret == 0) continue;
		else if (ret == -1) { client_goto_source = 1; goto client_thread_exit; }

		//Read data from global tsmpool and write it to client socket
		if(NMUX_DEBUG) fprintf(stderr, "client %p: sending...", param);
		ret = send(this_client->socket, pool_read_buffer + client_buffer_index, lpool->size - client_buffer_index, MSG_NOSIGNAL);
		if(NMUX_DEBUG) fprintf(stderr, "client sent.\n");
		if(ret == -1) 
		{
			switch(errno)
			{
				case EAGAIN: break;
				default: client_goto_source = 2; goto client_thread_exit;
			}
		}
		else client_buffer_index += ret;
	}

client_thread_exit:
	fprintf(stderr, "client %p: CS_THREAD_FINISHED, client_goto_source = %d, errno = %d", param, client_goto_source, errno);
	this_client->status = CS_THREAD_FINISHED;
	pthread_exit(NULL);
	return NULL;
}


int set_nonblocking(int fd)
{
	int flagtmp;
	if((flagtmp = fcntl(fd, F_GETFL))!=-1)
		if((flagtmp = fcntl(fd, F_SETFL, flagtmp|O_NONBLOCK))!=-1)
			return 0;
	return 1;
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
