#include <math.h>
#include "libcsdr.h"
#include "libcsdr_gpl.h"

typedef struct fastddc_s
{
	int pre_decimation;
	int post_decimation;
	int taps_length; 
	int taps_min_length;
	int overlap_length; //it is taps_length - 1
	int fft_size;
	int fft_inv_size;
	int input_size;
	int post_input_size;
	float pre_shift;
	int startbin; //for pre_shift
	int v; //step for pre_shift
	int offsetbin;
	float post_shift;
	int output_scrape;
	int scrap;
	shift_addition_data_t dsadata;
} fastddc_t;

int fastddc_init(fastddc_t* ddc, float transition_bw, int decimation, float shift_rate);
decimating_shift_addition_status_t fastddc_inv_cc(complexf* input, complexf* output, fastddc_t* ddc, FFT_PLAN_T* plan_inverse, complexf* taps_fft, decimating_shift_addition_status_t shift_stat);
void fastddc_print(fastddc_t* ddc, char* source);
void fft_swap_sides(complexf* io, int fft_size);
