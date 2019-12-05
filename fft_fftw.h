#pragma once

#ifdef USE_FFTW
//http://www.fftw.org/doc/Complex-One_002dDimensional-DFTs.html
//http://www.fftw.org/doc/Precision.html

#include <fftw3.h>
#define FFT_LIBRARY_USED "fftw3"

#define FFT_PLAN_T struct fft_plan_s
#define fft_malloc fftwf_malloc
#define fft_free fftwf_free

struct fft_plan_s 
{
	int size;
	void* input;
	void* output;
	fftwf_plan plan;
};

#include "libcsdr.h"

FFT_PLAN_T* make_fft_c2c(int size, complexf* input, complexf* output, int forward, int benchmark);
FFT_PLAN_T* make_fft_r2c(int size, float* input, complexf* output, int benchmark);
void fft_execute(FFT_PLAN_T* plan);
void fft_destroy(FFT_PLAN_T* plan);

/*
 * FFTW_MEASURE is inacceptably slow when there is no hardware cycle counter
 * available. Unfortunately, there's no way to detect this at compile- or
 * runtime.
 *
 * CSDR_DISABLE_FFTW_MEASURE can therefore be used to disable the use of
 * FFTW_MEASURE globally.
 *
 * additional information: http://www.fftw.org/fftw3_doc/Cycle-Counters.html
 *
 * https://github.com/simonyiszk/openwebrx/issues/139
 */
#ifdef CSDR_DISABLE_FFTW_MEASURE
#define CSDR_FFTW_MEASURE FFTW_ESTIMATE
#else
#define CSDR_FFTW_MEASURE FFTW_MEASURE
#endif

#endif
