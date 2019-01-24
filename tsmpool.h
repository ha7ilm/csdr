//tsmpool stands for Thread-Safe Memory Pool.

//It implements a big circular buffer that one thread writes into, and multiple threads read from.
//The reader threads have lower priority than the writer thread (they can be left behind if the don't read fast enough).

#include <vector>
#include <pthread.h>

#define TSM_DEBUG 0
#include <stdio.h>

using namespace std;

typedef struct tsmthread_s
{
	int read_index; //it always points to the next buffer to be read
} tsmthread_t;

class tsmpool
{
private:
	vector<tsmthread_t*> threads;
	vector<void*> buffers;
	int threads_cntr;
	pthread_mutex_t mutex;
	int ok; //tsmpool is expected to be included in C-style programs. 
			//	If something fails in the constructor, it will be seen here instead of a try{}catch{}
	int write_index; //it always points to the next buffer to be written
	int lowest_read_index; //unused
	int my_read_index; //it is used when tsmpool is used as a single writer - single reader circular buffer

public:
	const size_t size;
	const int num;
	int is_ok();
	tsmpool(size_t size, int num);
	void* get_write_buffer();
	tsmthread_t* register_thread();
	void remove_thread(tsmthread_t* thread);
	void* get_read_buffer(tsmthread_t* thread);
	int index_next(int index) { return (index+1==num)?0:index+1; }
	int index_before(int index) { return (index-1<0)?num-1:index-1; }
};
