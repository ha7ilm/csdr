/*
This software is part of libcsdr, a set of simple DSP routines for 
Software Defined Radio.

Copyright (c) 2014, Andras Retzler <randras@sdr.hu>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the copyright holder nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ANDRAS RETZLER BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once
#define MIN_M(x,y) (((x)>(y))?(y):(x))

/*
   _____                      _                                                 
  / ____|                    | |                                                
 | |     ___  _ __ ___  _ __ | | _____  __ 
 | |    / _ \| '_ ` _ \| '_ \| |/ _ \ \/ / 
 | |___| (_) | | | | | | |_) | |  __/>  <  
  \_____\___/|_| |_| |_| .__/|_|\___/_/\_\ 
                       | |                  
                       |_|                  
*/

typedef struct complexf_s { float i; float q; } complexf;

//apply to pointers:
#define iof(complexf_input_p,i) (*(((float*)complexf_input_p)+2*(i)))
#define qof(complexf_input_p,i) (*(((float*)complexf_input_p)+2*(i)+1))
#define absof(complexf_input_p,i) (sqrt(iof(complexf_input_p,i)*iof(complexf_input_p,i)+qof(complexf_input_p,i)*qof(complexf_input_p,i)))
#define argof(complexf_input_p,i) (atan2(qof(complexf_input_p,i),iof(complexf_input_p,i)))
#define cmult(cfo, cfi1, cfi2) iof(cfo,0)=iof(cfi1,0)*iof(cfi2,0)-qof(cfi1,0)*qof(cfi2,0);qof(cfo,0)=iof(cfi1,0)*qof(cfi2,0)+iof(cfi2,0)*qof(cfi1,0)
//(ai+aq*j)*(bi+bq*j)=ai*bi-aq*bq+(aq*bi+ai*bq)*j
#define cmultadd(cfo, cfi1, cfi2) { iof(cfo,0)+=iof(cfi1,0)*iof(cfi2,0)-qof(cfi1,0)*qof(cfi2,0);qof(cfo,0)+=iof(cfi1,0)*qof(cfi2,0)+iof(cfi2,0)*qof(cfi1,0); }
#define csetnull(cf) { iof(cf,0)=0.0; qof(cf,0)=0.0; }
#define e_powj(cf,w) { iof(cf,0)=cos(w); qof(cf,0)=sin(w); }
#define ccopy(dst,src) { iof(dst,0)=iof(src,0); qof(dst,0)=qof(src,0); }

//apply to values
#define iofv(complexf_input) (*((float*)&complexf_input))
#define qofv(complexf_input) (*(((float*)&complexf_input)+1))

//they dropped M_PI in C99, so we define it:
#define PI ((float)3.14159265358979323846)

//window
typedef enum window_s 
{
	WINDOW_BOXCAR, WINDOW_BLACKMAN, WINDOW_HAMMING
} window_t;

#define WINDOW_DEFAULT WINDOW_HAMMING

//FFT
//Note: these reference to things in this file (e.g. complexf):
#include "fft_fftw.h"
#include "fft_rpi.h"

// =================================================================================

//filter design
void firdes_lowpass_f(float *output, int length, float cutoff_rate, window_t window);
void firdes_bandpass_c(complexf *output, int length, float lowcut, float highcut, window_t window);
float firdes_wkernel_blackman(float input);
float firdes_wkernel_hamming(float input);
float firdes_wkernel_boxcar(float input);
window_t firdes_get_window_from_string(char* input);
char* firdes_get_string_from_window(window_t window);
int firdes_filter_len(float rolloff);

//demodulators
complexf fmdemod_quadri_cf(complexf* input, float* output, int input_size, float *temp, complexf last_sample);
complexf fmdemod_quadri_novect_cf(complexf* input, float* output, int input_size, complexf last_sample);
float fmdemod_atan_cf(complexf* input, float *output, int input_size, float last_phase);
void amdemod_cf(complexf* input, float *output, int input_size);
void amdemod_estimator_cf(complexf* input, float *output, int input_size, float alpha, float beta);
void limit_ff(float* input, float* output, int input_size, float max_amplitude);

//filters, decimators, resamplers, shift, etc.
int fir_decimate_cc(complexf *input, complexf *output, int input_size, int decimation, float *taps, int taps_length);
int deemphasis_nfm_ff (float* input, float* output, int input_size, int sample_rate);
float deemphasis_wfm_ff (float* input, float* output, int input_size, float tau, int sample_rate, float last_output);
float shift_math_cc(complexf *input, complexf* output, int input_size, float rate, float starting_phase);

typedef struct dcblock_preserve_s
{
	float last_input;
	float last_output;
} dcblock_preserve_t;
dcblock_preserve_t dcblock_ff(float* input, float* output, int input_size, float a, dcblock_preserve_t preserved);
float fastdcblock_ff(float* input, float* output, int input_size, float last_dc_level);

typedef struct fastagc_ff_s
{
	float* buffer_1; 
	float* buffer_2;
	float* buffer_input; //it is the actual input buffer to fill
	float peak_1;
	float peak_2;
	int input_size;
	float reference;
	float last_gain;
} fastagc_ff_t;

void fastagc_ff(fastagc_ff_t* input, float* output);

typedef struct rational_resampler_ff_s
{
	int input_processed;
	int output_size;
	int last_taps_delay;
} rational_resampler_ff_t;

rational_resampler_ff_t rational_resampler_ff(float *input, float *output, int input_size, int interpolation, int decimation, float *taps, int taps_length, int last_taps_delay);
void rational_resampler_get_lowpass_f(float* output, int output_size, int interpolation, int decimation, window_t window);

void apply_window_c(complexf* input, complexf* output, int size, window_t window);
void apply_window_f(float* input, float* output, int size, window_t window);
void logpower_cf(complexf* input, float* output, int size, float add_db);

typedef struct fractional_decimator_ff_s
{
	float remain;
	int input_processed;
	int output_size;
} fractional_decimator_ff_t;
fractional_decimator_ff_t fractional_decimator_ff(float* input, float* output, int input_size, float rate, float *taps, int taps_length, fractional_decimator_ff_t d);

int log2n(int x);
int next_pow2(int x);
void apply_fir_fft_cc(FFT_PLAN_T* plan, FFT_PLAN_T* plan_inverse, complexf* taps_fft, complexf* last_overlap, int overlap_size);
void gain_ff(float* input, float* output, int input_size, float gain);

void convert_u8_f(unsigned char* input, float* output, int input_size);
void convert_f_u8(float* input, unsigned char* output, int input_size);
void convert_f_i16(float* input, short* output, int input_size);
void convert_i16_f(short* input, float* output, int input_size);




