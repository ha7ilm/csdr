#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <signal.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <vector>
#include <limits.h>

#define SOFTWARE_NAME "ddcd"
#define MSG_START SOFTWARE_NAME ": "

typedef enum ddc_method_e
{
	M_TD,
	M_FASTDDC
} ddc_method_t;

typedef enum client_status_e
{
	CS_CREATED,
	CS_THREAD_RUNNING,
	CS_THREAD_FINISHED
} client_status_t;


typedef struct client_s
{
	struct sockaddr_in addr;
	int socket;
	int error; //set to non-zero on error (data transfer failed)
	pthread_t thread;
	client_status_t status;

} client_t;

typedef enum command_type_e
{
	CT_SHIFT,
	CT_BYPASS
} command_type_t;


typedef struct command_s
{
	command_type_t type;
	float float_param;
} command_t;

void print_exit(const char* why);
void error_exit(const char* why);
void maxfd(int* maxfd, int fd);
