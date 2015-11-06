#pragma once

#include <signal.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <iostream>
#include <vector>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/wait.h>

typedef struct client_s
{
	struct sockaddr_in addr;
	int socket;
	pid_t pid;
	int pipefd[2];
	int error;
} client_t;


void client();
void error_exit(const char* why);
void print_client(client_t* client, const char* what);
int proc_exists(pid_t pid);
