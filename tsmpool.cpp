tsmpool::tsmpool(size_t size, int num) 
{
	this->threads_cntr = 0;
	this->num = num;
	this->size = size;
	this->ok = 1;
	this->lowest_read_index = -1;
    if (pthread_mutex_init(&this->mutex, NULL) != 0) this->ok=0;
}

size_t tsmpool::get_size() { return this->size; }

void* tsmpool::get_write_buffer() 
{  
	if(write_index==index_before(lowest_read_index)) return NULL;
	void* to_return = buffers[write_index];
	write_index=index_next(write_index);
}

tsmthread_t* tsmpool::register_thread()
{
	if(!ok) return -1;
	pthread_mutex_lock(&this->mutex);
	tsmthread_t* thread = new tsmthread_t;
	thread->read_index = write_index;
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
	if(thread->read_index==write_index) return NULL;
	void* to_return = buffers[thread->read_index];
	thread->read_index=index_next(thread->read_index);
}

void* tsmpool::set_read_index_distance(tsmthread_t* thread, int distance)
{
}
