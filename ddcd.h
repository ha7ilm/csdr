#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#define SOFTWARE_NAME "ddcd"
#define MSG_START SOFTWARE_NAME ": "

typedef enum ddc_method_e 
{
	M_TD,
	M_FASTDDC
} ddc_method_t;

void print_exit(const char* why);

typedef struct client_s
{
	struct sockaddr_in addr;
	int socket;
	int error;
	pthread_t thread;
} client_t;

