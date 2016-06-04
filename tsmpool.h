//tsmpool stands for Thread-Safe Memory Pool.

//It implements a big circular buffer that one thread writes into, and multiple threads read from.
//The reader threads have lower priority than the writer thread (they can be left behind if the don't read fast enough).

typedef struct tsmthread_s
{
	int read_index; //it always points to the next buffer to be read
} tsmthread_t;

class tsmpool
{
private:
	size_t size;
	int num;
	vector<tsmthread_t*> threads;
	vector<void*> buffers;
	int threads_cntr;
	pthread_mutex_t mutex;
	int ok;
	int write_index; //it always points to the next buffer to be written
	int lowest_read_index; //unused

public:
	size_t get_size();
	tsmpool(size_t size, int num);
	void* get_write_buffer();
	int register_thread();
	void* get_read_buffer(int thread_id);
	int index_next(int index) { return (index+1==size)?0:index; }
	int index_before(int index) { return (index-1<0)?size-1:index; }
}
