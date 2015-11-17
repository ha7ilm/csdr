#include <math.h>
#include "libcsdr.h"
#include "libcsdr_gpl.h"

typedef struct fastddc_s
{
	int pre_decimation;
	int post_decimation;
	int taps_length; 
	int taps_real_length;
	int overlap_length; //it is taps_length - 1
	int fft_size;
	int fft_inv_size;
	int input_size;
	int output_size;
	float pre_shift;
	int startbin; //for pre_shift
	int v; //step for pre_shift
	int offsetbin;
	float post_shift;
	int output_scrape;
	int scrape;
} fastddc_t;

int fastddc_init(fastddc_t* ddc, float transition_bw, int decimation, float shift_rate);
decimating_shift_addition_status_t fastddc_inv_cc(complexf* input, complexf* output, fastddc_t* ddc, FFT_PLAN_T* plan_inverse, complexf* taps_fft, decimating_shift_addition_status_t shift_stat);
void fastddc_print(fastddc_t* ddc);
