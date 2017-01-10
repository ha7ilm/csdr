#include "tsmpool.h"

tsmpool::tsmpool(size_t size, int num)
	size(size), 
	num(num) //number of buffers of (size) to alloc
{
	this->threads_cntr = 0;
	this->ok = 1;
	this->lowest_read_index = -1;
	this->write_index = 0;
	this->my_read_index = 0;
    if (pthread_mutex_init(&this->mutex, NULL) != 0) { this->ok = 0; return; }
	for(int i=0; i<num; i++) 
	{
		void* newptr = (void*)new char[size];
		if(!newptr) { this->ok = 0; return; }
		buffers.push_back(newptr);
	}
}

int tsmpool::is_ok() { return this->ok; }

void* tsmpool::get_write_buffer()
{
	//if(write_index==index_before(lowest_read_index)) return NULL;
	void* to_return = buffers[write_index];
	write_index=index_next(write_index);
}

tsmthread_t* tsmpool::register_thread()
{
	if(!ok) return NULL;
	pthread_mutex_lock(&this->mutex);
	tsmthread_t* thread = new tsmthread_t;
	thread->read_index = index_before(write_index);
	threads.push_back(thread);
	pthread_mutex_unlock(&this->mutex);
	return thread;
}

int tsmpool::remove_thread(tsmthread_t* thread)
{
	pthread_mutex_lock(&this->mutex);
	for(int i=0;i<threads.size();i++)
		if(threads[i] == thread)
		{
			delete threads[i];
			threads.erase(i);
			break;
		}
	pthread_mutex_unlock(&this->mutex);
}

void* tsmpool::get_read_buffer(tsmthread_t* thread)
{
	int* actual_read_index = (thread==NULL) ? &my_read_index : &thread->read_index;
	if(*actual_read_index==index_before(write_index)) return NULL;
	void* to_return = buffers[*actual_read_index];
	*actual_read_index=index_next(*actual_read_index);
}
