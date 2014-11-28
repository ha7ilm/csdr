#ifdef USE_FFTW

#include "fft_fftw.h"
#include <stdlib.h>

FFT_PLAN_T* make_fft_c2c(int size, complexf* input, complexf* output, int forward, int benchmark)
{
	FFT_PLAN_T* plan=(FFT_PLAN_T*)malloc(sizeof(FFT_PLAN_T));
	plan->plan = fftwf_plan_dft_1d(size, (fftwf_complex*)input, (fftwf_complex*)output, (forward)?FFTW_FORWARD:FFTW_BACKWARD, (benchmark)?FFTW_MEASURE:FFTW_ESTIMATE);
	plan->size=size;
	plan->input=(void*)input;
	plan->output=(void*)output;
	return plan;
}

FFT_PLAN_T* make_fft_r2c(int size, float* input, complexf* output, int benchmark) //always forward DFT
{
	FFT_PLAN_T* plan=(FFT_PLAN_T*)malloc(sizeof(FFT_PLAN_T));
	plan->plan = fftwf_plan_dft_r2c_1d(size, input, (fftwf_complex*)output, (benchmark)?FFTW_MEASURE:FFTW_ESTIMATE);
	plan->size=size;
	plan->input=(void*)input;
	plan->output=(void*)output;
	return plan;
}

FFT_PLAN_T* make_fft_c2r(int size, complexf* input, float* output, int benchmark) //always backward DFT
{
	FFT_PLAN_T* plan=(FFT_PLAN_T*)malloc(sizeof(FFT_PLAN_T));
	plan->plan = fftwf_plan_dft_c2r_1d(size, (fftwf_complex*)input, output, (benchmark)?FFTW_MEASURE:FFTW_ESTIMATE);
	plan->size=size;
	plan->input=(void*)input;
	plan->output=(void*)output;
	return plan;
}

void fft_execute(FFT_PLAN_T* plan)
{
	fftwf_execute(plan->plan);
}

void fft_destroy(FFT_PLAN_T* plan)
{
	fftwf_destroy_plan(plan->plan);
	free(plan);
}

#endif
