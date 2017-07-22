#include "wspr.h"

// used by qsort
// NB: assumes first element in struct is float sort field
int qsort_floatcomp(const void *elem1, const void *elem2)
{
	const float f1 = *(const float *) elem1, f2 = *(const float *) elem2;
    if (f1 < f2)
        return -1;
    return f1 > f2;
}

static bool init = false;
static u4_t epoch_sec;

static void set_epoch()
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	epoch_sec = ts.tv_sec;
	
	init = true;
}

// overflows 136 years after timer epoch
u4_t timer_sec()
{
	struct timespec ts;

	if (!init) set_epoch();
	clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec - epoch_sec;
}

// overflows 49.7 days after timer epoch
u4_t timer_ms()
{
	struct timespec ts;

	if (!init) set_epoch();
	clock_gettime(CLOCK_MONOTONIC, &ts);
	int msec = ts.tv_nsec/1000000;
	assert(msec >= 0 && msec < 1000);
    return (ts.tv_sec - epoch_sec)*1000 + msec;
}
