#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <signal.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <errno.h>
#include <fcntl.h>
#include <poll.h>
#include <arpa/inet.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include "tsmpool.h"

#define MSG_START "nmux: "
#define NMUX_DEBUG 0

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
    tsmthread_t* tsmthread;
	client_status_t status;
    //the following members are there to give access to some global variables inside the thread:
    tsmpool* lpool; 
    int sleeping;
} client_t;

void print_exit(const char* why);
void sig_handler(int signo);
void* client_thread (void* param);
void error_exit(const char* why);
void maxfd(int* maxfd, int fd);
int set_nonblocking(int fd);
