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
float transition_bw = 0.05;
int bufsize = 1024;
int bufsizeall;
int pipe_max_size;
int in_client = 0;
char ddc_method_str[100] = "td";
ddc_method_t ddc_method;
pid_t main_dsp_proc;

int input_fd = STDIN_FILENO; //can be stdin, or the stdout of main_subprocess
pid_t main_subprocess_pid = 0;
pid_t main_subprocess_pgrp = 0;
pid_t client_subprocess_pid = 0;
pid_t client_subprocess_pgrp = 0;

char* buf;

int set_nonblocking(int fd)
{
	int flagtmp;
	if((flagtmp = fcntl(fd, F_GETFL))!=-1)
		if((flagtmp = fcntl(fd, F_SETFL, flagtmp|O_NONBLOCK))!=-1)
			return 0;
	return 1;
}

int proc_exists(pid_t pid)
{
	if(pid==0 || pid==1) return 1;
	return kill(pid, 0) != -1;
}

void sig_handler(int signo)
{	
	int tmpstat;
	if(signo==SIGPIPE) 
	{
		fprintf(stderr,MSG_START "SIGPIPE received.\n");
		return;
	}	
	if(signo==SIGCHLD)
		if(  main_subprocess_pid  && signo==SIGCHLD && (waitpid(main_subprocess_pid, &tmpstat, WNOHANG), 1) && !proc_exists(main_subprocess_pid)  )
		{
			fprintf(stderr,MSG_START "main_subprocess_pid exited! Exiting...\n");
		}
	else return;
	//if(pgrp!=1 && pgrp!=0) //I just want to make sure that we cannot kill init or sched
	//	killpg(pgrp, signo);
	if( !in_client && main_subprocess_pid ) killpg2(main_subprocess_pgrp);
	if( in_client && client_subprocess_pid ) killpg2(client_subprocess_pgrp);
	fprintf(stderr, MSG_START "signal %d caught in %s, exiting ddcd...\n", signo, (in_client)?"client":"main");
	fflush(stderr);
	exit(0);
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

	//set signals
	struct sigaction sa;
	memset(&sa, 0, sizeof(sa));
	sa.sa_handler = sig_handler;
	sigaction(SIGTERM, &sa, NULL);
	sigaction(SIGKILL, &sa, NULL);
	sigaction(SIGQUIT, &sa, NULL);
	sigaction(SIGINT, &sa, NULL);
	sigaction(SIGHUP, &sa, NULL);
	sigaction(SIGCHLD, &sa, NULL);
	sigaction(SIGPIPE, &sa, NULL);
	prctl(PR_SET_PDEATHSIG, SIGHUP); //get a signal when parent exits

	struct sockaddr_in addr_host;
    int listen_socket;
	std::vector<client_t*> clients;
	clients.reserve(100);
    listen_socket=socket(AF_INET,SOCK_STREAM,0);

	int sockopt = 1;
	if( setsockopt(listen_socket, SOL_SOCKET, SO_REUSEADDR, (char *)&sockopt, sizeof(sockopt)) == -1 )
		error_exit(MSG_START "cannot set SO_REUSEADDR");  //the best description on SO_REUSEADDR ever: http://stackoverflow.com/a/14388707/3182453

    memset(&addr_host,'0',sizeof(addr_host));
    addr_host.sin_family=AF_INET;
    addr_host.sin_port=htons(host_port);
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

	bufsizeall = bufsize*sizeof(char);
	buf = (char*)malloc(bufsizeall);

	FILE* tempfile = fopen("/proc/sys/fs/pipe-max-size","r");
	if(!tempfile)
	{
		perror(MSG_START "cannot read /proc/sys/fs/pipe-max-size");
	}
	else
	{
		char pipe_max_size_str[100];
		int tfread = fread(pipe_max_size_str, 1, 100, tempfile);
		pipe_max_size_str[tfread]='\0';
		pipe_max_size = atoi(pipe_max_size_str);
		//fprintf(stderr, MSG_START "note: pipe_max_size = %d\n", pipe_max_size);
		//if(pipe_max_size>4096 && fcntl(STDIN_FILENO, F_SETPIPE_SZ, pipe_max_size)==-1)
		//	perror("failed to fcntl(STDIN_FILENO, F_SETPIPE_SZ, ...)");
	}

	//We'll see if it is a good idea:
	//setpgrp();
	//pgrp = getpgrp();
	//It is not, because we can't catch Ctrl+C (SIGINT), as it is sent to a process group...

	//Start DSP subprocess from the main process if required
	char main_subprocess_cmd_buf[500];


	int pipe_m2s_ctl[2];	//main to subprocess :: control channel
	int pipe_s2m[2];		//subprocess to main

	if(pipe(pipe_m2s_ctl)) error_exit(MSG_START "couldn't create pipe_m2s_ctl");
	if(pipe(pipe_s2m)) error_exit(MSG_START "couldn't create pipe_s2m");

	if(decimation!=1)
	{
		switch(ddc_method)
		{
		case M_TD:
			break;
		case M_FASTDDC:
			sprintf(main_subprocess_cmd_buf, subprocess_args_fastddc_1, decimation, transition_bw);
			fprintf(stderr, MSG_START "starting main_subprocess_cmd: %s\n", main_subprocess_cmd_buf);
			if(!(main_subprocess_pid = run_subprocess( main_subprocess_cmd_buf, 0, pipe_s2m, &main_subprocess_pgrp )))
				print_exit(MSG_START "couldn't start main_subprocess_cmd!\n");
			close(STDIN_FILENO); // redirect stdin to the stdin of the subprocess 
			break;
		}
	}

	int highfd = 0; 
	FD_ZERO(&select_fds);
	FD_SET(listen_socket, &select_fds);
	maxfd(&highfd, listen_socket);
	if(main_subprocess_pid) input_fd = pipe_s2m[0]; //else STDIN_FILENO
	FD_SET(input_fd, &select_fds);
	maxfd(&highfd, input_fd);

	//Set stdin and listen_socket to non-blocking 
	if(set_nonblocking(input_fd) || set_nonblocking(listen_socket)) //don't do it before subprocess fork!
		error_exit(MSG_START "cannot set_nonblocking()");

	for(;;)
	{
		//Let's wait until there is any new data to read, or any new connection!
		select(highfd, &select_fds, NULL, NULL, NULL);

		//Is there a new client connection?
		if( (new_socket = accept(listen_socket, (struct sockaddr*)&addr_cli, &addr_cli_len)) != -1)
		{ 
			this_client = new client_t;
			this_client->error = 0;
			memcpy(&this_client->addr, &addr_cli, sizeof(this_client->addr));
			this_client->socket = new_socket;
			if(pipe(this_client->pipefd) == -1)
			{ 
				perror(MSG_START "cannot open new pipe() for the client");
				continue;
			}
			if(fcntl(this_client->pipefd[1], F_SETPIPE_SZ, pipe_max_size) == -1)
				perror("failed to F_SETPIPE_SZ for the client pipe");
			if(this_client->pid = fork())
			{
				//We're the parent
				set_nonblocking(this_client->pipefd[1]);
				clients.push_back(this_client);
				fprintf(stderr, MSG_START "client pid: %d\n", this_client->pid);
			}
			else
			{
				//We're the client
				client();
				return 1;
			}
		}

		int retval = read(input_fd, buf, bufsizeall);
		if(retval==0)
		{
			//end of input stream, close clients and exit
		}
		else if(retval != -1)
		{
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
	}

	return 0; 
}

pid_t run_subprocess(char* cmd, int* pipe_in, int* pipe_out, pid_t* pgrp)
{
	/*char sem_name[101];
	snprintf(sem_name,100,"ddcd_sem_%d",getpid());
	sem_t mysem;
	if(sem_init(&mysem, 1, 1)==-1) error_exit("failed to sem_init() in run_subprocess()");
	fprintf(stderr, "sem_waiting\n");
	if(sem_wait(&mysem)==-1) error_exit("the first sem_wait() failed in run_subprocess()");
	fprintf(stderr, "sem_waited\n");
	*/
	int syncpipe[2];
	if(pipe(syncpipe)==-1) error_exit("failed to create pipe()");
	pid_t pid = fork();

	if(pid < 0) return 0; //fork failed
	if(pid == 0)
	{
		setpgrp();
		write(syncpipe[1], " ", 1);
		//if(sem_post(&mysem)==-1) error_exit("failed to sem_post() in run_subprocess()");
		//We're the subprocess
		//fprintf(stderr, "run_subprocess :: execl\n");
		//if(fcntl(pipe_in[1], F_SETPIPE_SZ, pipe_max_size) == -1) perror("Failed to F_SETPIPE_SZ in run_subprocess()");
		if(pipe_in)
		{
			close(pipe_in[1]);
			dup2(pipe_in[0], STDIN_FILENO);
		}
		if(pipe_out)
		{
			close(pipe_out[0]);
			dup2(pipe_out[1], STDOUT_FILENO);
		}
		execl("/bin/bash","bash","-c",cmd, (char*)0);
		error_exit(MSG_START "run_subprocess failed to execute command");
	}
	else
	{ 
		//if(sem_wait(&mysem)==-1) error_exit("the second sem_wait() failed in run_subprocess()");
		int synctemp;
		read(syncpipe[0], &synctemp, 1);
		*pgrp = getpgid(pid);
		fprintf(stderr, MSG_START "run_subprocess pgid returned = %d\n", *pgrp);
		return pid;
	}
}

void print_client(client_t* client, const char* what)
{
	fprintf(stderr,MSG_START "(client %s:%d) %s\n", inet_ntoa(client->addr.sin_addr), client->addr.sin_port, what);
}

#define CTL_BUFSIZE 1024

int read_socket_ctl(int fd, char* output, int max_size)
{
	//fprintf(stderr, "doing read_socket_ctl %d\n", fd);
	//if(!fd) return 0;
	static char buffer[CTL_BUFSIZE];
	static int buffer_index=0;
	if(buffer_index==CTL_BUFSIZE) buffer_index=0;
	int bytes_read=recv(fd,buffer+buffer_index,(CTL_BUFSIZE-buffer_index)*sizeof(char), MSG_DONTWAIT);
	if(bytes_read<=0) return 0;
	//fprintf(stderr, "recv %d\n", bytes_read);
	
	int prev_newline_at=0;
	int last_newline_at=0;
	for(int i=0;i<buffer_index+bytes_read;i++) 
	{
		if(buffer[i]=='\n') 
		{
			prev_newline_at=last_newline_at;
			last_newline_at=i+1;
		}
	}
	if(last_newline_at)
	{
		int oi=0;
		for(int i=prev_newline_at;buffer[i]!='\n'&&oi<max_size;i++) output[oi++]=buffer[i]; //copy to output buffer
		output[oi++]='\0';
		memmove(buffer,buffer+last_newline_at,buffer_index+bytes_read-last_newline_at);
		buffer_index=bytes_read-last_newline_at;
		return 1;
	}
	else
	{
		buffer_index+=bytes_read;
	 	return 0;
	}
}

int ctl_get_arg(char* input, const char* cmd, const char* format, ...)
{
	int retval=0;
	int cmdlen=strlen(cmd);
	if(input[cmdlen]=='=')
	{
		//fprintf(stderr, "cga found=\n");
		if(input[cmdlen]=0, !strcmp(input,cmd))
		{
			//fprintf(stderr, "cga foundokay\n");
			va_list vl;
			va_start(vl,format);
			retval=vsscanf(input+cmdlen+1,format,vl);
			va_end(vl);
		}
		input[cmdlen]='=';
	}
	//fprintf(stderr, "cga retval %d\n", retval);
	return retval;
}

void client()
{
	in_client=1;
	print_client(this_client, "client process forked.");
	
	char client_subprocess_cmd_buf[500];
	int input_fd = this_client->pipefd[0];
	int pipe_ctl[2], pipe_stdout[2];

	prctl(PR_SET_PDEATHSIG, SIGHUP); //get a signal when parent exits
	
	if(decimation!=1)
	{
		
		if(pipe(pipe_ctl)==-1) error_exit(MSG_START "cannot open new pipe() for the client subprocess");
		if(pipe(pipe_stdout)==-1) error_exit(MSG_START "cannot open new pipe() for the client subprocess");
		switch(ddc_method)
		{
		case M_TD:
			sprintf(client_subprocess_cmd_buf, subprocess_cmd_td, pipe_ctl[0], decimation, transition_bw);
			break;
		case M_FASTDDC:
			sprintf(client_subprocess_cmd_buf, subprocess_args_fastddc_2, pipe_ctl[0], decimation, transition_bw);			
			break;
		}

		if(!(client_subprocess_pid = run_subprocess( client_subprocess_cmd_buf, this_client->pipefd, pipe_stdout, &client_subprocess_pgrp))) 
			print_exit(MSG_START "couldn't start client_subprocess_cmd!\n");
		fprintf(stderr, MSG_START "starting client_subprocess_cmd: %s\n", client_subprocess_cmd_buf);
		input_fd = pipe_stdout[0]; //we don't have to set it nonblocking
		fprintf(stderr, MSG_START "pipe_stdout[0] = %d\n", pipe_stdout[0]);
		write(pipe_ctl[1], "0.0\n", 4);
	}
	char recv_cmd[CTL_BUFSIZE];
	char temps[CTL_BUFSIZE*2];
	int tempi;
	float tempf;

	for(;;)
	{
		while(read_socket_ctl(this_client->socket, recv_cmd, CTL_BUFSIZE)) 
		{
			sprintf(temps, "read_socket_ctl: %s", recv_cmd);
			print_client(this_client, temps);
			if(ctl_get_arg(recv_cmd, "bypass", "%d", &tempi))
			{
				if(tempi==1 && client_subprocess_pid)
				{
					//print_client(this_client, "suspending client_subprocess_pgrp...\n");
					//fprintf(stderr, "client_subprocess_pgrp = %d\n", client_subprocess_pgrp);
					//killpg(client_subprocess_pgrp, SIGTSTP);
					//while(proc_exists(client_subprocess_pid)) usleep(10000);
					//print_client(this_client, "done killing client_subprocess_pid.\n");
					input_fd=this_client->pipefd[0]; //by doing this, we don't read from pipe_stdout[0] anymore, so that csdr stops doing anything, and also doesn't read anymore from the input: we get the whole I/Q stream!
				}
				if(tempi==0 && client_subprocess_pid)
				{
					input_fd=pipe_stdout[0];
				}
				
			}
			if(ctl_get_arg(recv_cmd, "shift", "%g", &tempf))
			{			
				tempi=sprintf(temps, "%g\n", tempf);
				write(pipe_ctl[1], temps, tempi);
				fsync(pipe_ctl[1]);
			}
		}
		int nread = read(input_fd,buf,bufsizeall);
		if(nread<=0) continue;
		if(send(this_client->socket,buf,nread,0)==-1)
		{
			print_client(this_client, "client process is exiting.\n");
			if(client_subprocess_pid) killpg2(client_subprocess_pgrp);
			exit(0);
		}
	}	
}

void killpg2(pid_t pgrp)
{
	//fprintf(stderr, MSG_START "killpg2: %d\n", pgrp);
	if(pgrp!=1 && pgrp!=0) killpg(pgrp, SIGTERM);
}

void error_exit(const char* why)
{
	perror(why);
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
