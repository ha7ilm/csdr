#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <signal.h>
#include <arpa/inet.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <vector>

#define SOFTWARE_NAME "ddcd"
#define MSG_START SOFTWARE_NAME ": "

typedef enum ddc_method_e 
{
	M_TD,
	M_FASTDDC
} ddc_method_t;

typedef struct client_s
{
	struct sockaddr_in addr;
	int id;
	int socket;
	int error;
	pthread_t thread;
} client_t;

void print_exit(const char* why);
void error_exit(const char* why);
void maxfd(int* maxfd, int fd);


