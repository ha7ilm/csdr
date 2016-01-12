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
#include <sys/prctl.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <semaphore.h>

typedef struct client_s
{
	struct sockaddr_in addr;
	int socket;
	pid_t pid;
	int pipefd[2];
	int error;
	pid_t dsp_proc;
} client_t;


void client();
void error_exit(const char* why);
void print_exit(const char* why);
void print_client(client_t* client, const char* what);
int proc_exists(pid_t pid);
pid_t run_subprocess(char* cmd, int* pipe_in, int* pipe_out, pid_t* pgrp);
void maxfd(int* maxfd, int fd);
void sig_handler(int signo);
void killpg2(pid_t pgrp);
int ctl_get_arg(char* input, const char* cmd, const char* format, ...);

typedef enum ddc_method_e 
{
	M_TD,
	M_FASTDDC
} ddc_method_t;

const char subprocess_cmd_td[] = "csdr "
#ifdef NEON_OPTS
	"shift_addfast_cc"
#else
	"shift_unroll_cc"
#endif
	" --fd %d | csdr fir_decimate_cc %d %g";

const char subprocess_args_fastddc_1[] = "csdr fastddc_fwd_cc %d %g";
//const char subprocess_args_fastddc_1[] = "csdr through %d %g";
const char subprocess_args_fastddc_2[] = "csdr fastddc_inv_cc --fd %d %d %g";
//const char subprocess_args_fastddc_2[] = "csdr convert_u8_f %d %d %g";
