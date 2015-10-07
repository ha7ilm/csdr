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

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include "libcsdr.h"
#include "predefined.h"
#include <assert.h>

/*
           _           _                   __                  _   _                 
          (_)         | |                 / _|                | | (_)                
 __      ___ _ __   __| | _____      __  | |_ _   _ _ __   ___| |_ _  ___  _ __  ___ 
 \ \ /\ / / | '_ \ / _` |/ _ \ \ /\ / /  |  _| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
  \ V  V /| | | | | (_| | (_) \ V  V /   | | | |_| | | | | (__| |_| | (_) | | | \__ \
   \_/\_/ |_|_| |_|\__,_|\___/ \_/\_/    |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
                                                                                     
                                                                                                                                                         
*/

#define MFIRDES_GWS(NAME) \
	if(!strcmp( #NAME , input )) return WINDOW_ ## NAME;

window_t firdes_get_window_from_string(char* input)
{
	MFIRDES_GWS(BOXCAR);
	MFIRDES_GWS(BLACKMAN);
	MFIRDES_GWS(HAMMING);
	return WINDOW_DEFAULT;
}

#define MFIRDES_GSW(NAME) \
	if(window == WINDOW_ ## NAME) return #NAME;

char* firdes_get_string_from_window(window_t window)
{
	MFIRDES_GSW(BOXCAR);
	MFIRDES_GSW(BLACKMAN);
	MFIRDES_GSW(HAMMING);
	return "INVALID";
}

float firdes_wkernel_blackman(float rate)
{
	//Explanation at Chapter 16 of dspguide.com, page 2
	//Blackman window has better stopband attentuation and passband ripple than Hamming, but it has slower rolloff.
	rate=0.5+rate/2;
	return 0.42-0.5*cos(2*PI*rate)+0.08*cos(4*PI*rate);
}

float firdes_wkernel_hamming(float rate)
{
	//Explanation at Chapter 16 of dspguide.com, page 2
	//Hamming window has worse stopband attentuation and passband ripple than Blackman, but it has faster rolloff.
	rate=0.5+rate/2;
	return 0.54-0.46*cos(2*PI*rate);
}


float firdes_wkernel_boxcar(float rate)
{	//"Dummy" window kernel, do not use; an unwindowed FIR filter may have bad frequency response
	return 1.0;
}

float (*firdes_get_window_kernel(window_t window))(float)
{	
	if(window==WINDOW_HAMMING) return firdes_wkernel_hamming;
	else if(window==WINDOW_BLACKMAN) return firdes_wkernel_blackman;
	else if(window==WINDOW_BOXCAR) return firdes_wkernel_boxcar;
	else return firdes_get_window_kernel(WINDOW_DEFAULT);
}

/*
  ______ _____ _____      __ _ _ _                   _           _             
 |  ____|_   _|  __ \    / _(_) | |                 | |         (_)            
 | |__    | | | |__) |  | |_ _| | |_ ___ _ __     __| | ___  ___ _  __ _ _ __  
 |  __|   | | |  _  /   |  _| | | __/ _ \ '__|   / _` |/ _ \/ __| |/ _` | '_ \ 
 | |     _| |_| | \ \   | | | | | ||  __/ |     | (_| |  __/\__ \ | (_| | | | |
 |_|    |_____|_|  \_\  |_| |_|_|\__\___|_|      \__,_|\___||___/_|\__, |_| |_|
                                                                    __/ |      
                                                                   |___/       
*/

void firdes_lowpass_f(float *output, int length, float cutoff_rate, window_t window)
{	//Generates symmetric windowed sinc FIR filter real taps
	//	length should be odd
	//	cutoff_rate is (cutoff frequency/sampling frequency)
	//Explanation at Chapter 16 of dspguide.com
	int middle=length/2;
	float temp;
	float (*window_function)(float)  = firdes_get_window_kernel(window);
	output[middle]=2*PI*cutoff_rate*window_function(0);
	for(int i=1; i<=middle; i++) //@@firdes_lowpass_f: calculate taps
	{
		output[middle-i]=output[middle+i]=(sin(2*PI*cutoff_rate*i)/i)*window_function((float)i/middle);
		//printf("%g %d %d %d %d | %g\n",output[middle-i],i,middle,middle+i,middle-i,sin(2*PI*cutoff_rate*i));
	}
	
	//Normalize filter kernel
	float sum=0;
	for(int i=0;i<length;i++) //@firdes_lowpass_f: normalize pass 1
	{
		sum+=output[i];
	}
	for(int i=0;i<length;i++) //@firdes_lowpass_f: normalize pass 2
	{
		output[i]/=sum;
	}
}

void firdes_bandpass_c(complexf *output, int length, float lowcut, float highcut, window_t window)
{
	//To generate a complex filter:
	//	1. we generate a real lowpass filter with a bandwidth of highcut-lowcut
	//	2. we shift the filter taps spectrally by multiplying with e^(j*w), so we get complex taps
	//(tnx HA5FT)
	float* realtaps = (float*)malloc(sizeof(float)*length);
	
	firdes_lowpass_f(realtaps, length, (highcut-lowcut)/2, window);
	float filter_center=(highcut+lowcut)/2;
	
	float phase=0, sinval, cosval;
	for(int i=0; i<length; i++) //@@firdes_bandpass_c
	{
		cosval=cos(phase);
		sinval=sin(phase);		
		phase+=2*PI*filter_center;
		while(phase>2*PI) phase-=2*PI; //@@firdes_bandpass_c
		while(phase<0) phase+=2*PI; 
		iof(output,i)=cosval*realtaps[i];
		qof(output,i)=sinval*realtaps[i];
		//output[i] := realtaps[i] * e^j*w
	}
}

int firdes_filter_len(float transition_bw)
{
	int result=4.0/transition_bw;
	if (result%2==0) result++; //number of symmetric FIR filter taps should be odd
	return result;
}

/*
  _____   _____ _____      __                  _   _                 
 |  __ \ / ____|  __ \    / _|                | | (_)                
 | |  | | (___ | |__) |  | |_ _   _ _ __   ___| |_ _  ___  _ __  ___ 
 | |  | |\___ \|  ___/   |  _| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
 | |__| |____) | |       | | | |_| | | | | (__| |_| | (_) | | | \__ \
 |_____/|_____/|_|       |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___/

*/                                                   

float shift_math_cc(complexf *input, complexf* output, int input_size, float rate, float starting_phase)
{
	rate*=2;
	//Shifts the complex spectrum. Basically a complex mixer. This version uses cmath.
	float phase=starting_phase;
	float phase_increment=rate*PI;
	float cosval, sinval;
	for(int i=0;i<input_size; i++) //@shift_math_cc
	{
		cosval=cos(phase);
		sinval=sin(phase);
		//we multiply two complex numbers.
		//how? enter this to maxima (software) for explanation:
		//   (a+b*%i)*(c+d*%i), rectform;
		iof(output,i)=cosval*iof(input,i)-sinval*qof(input,i);
		qof(output,i)=sinval*iof(input,i)+cosval*qof(input,i);
		phase+=phase_increment;
		while(phase>2*PI) phase-=2*PI; //@shift_math_cc: normalize phase
		while(phase<0) phase+=2*PI;
	}
	return phase;
}



shift_table_data_t shift_table_init(int table_size)
{
	//RTODO
	shift_table_data_t output;
	output.table=(float*)malloc(sizeof(float)*table_size);
	output.table_size=table_size;
	for(int i=0;i<table_size;i++)
	{
		output.table[i]=sin(((float)i/table_size)*(PI/2));
	}
	return output;
}

void shift_table_deinit(shift_table_data_t table_data)
{
	free(table_data.table);
}

float shift_table_cc(complexf* input, complexf* output, int input_size, float rate, shift_table_data_t table_data, float starting_phase)
{
	//RTODO
	rate*=2;
	//Shifts the complex spectrum. Basically a complex mixer. This version uses a pre-built sine table.
	float phase=starting_phase;
	float phase_increment=rate*PI;
	float cosval, sinval;
	for(int i=0;i<input_size; i++) //@shift_math_cc
	{
		int sin_index, cos_index, temp_index, sin_sign, cos_sign;
		//float vphase=fmodf(phase,PI/2); //between 0 and 90deg
		int quadrant=phase/(PI/2); //between 0 and 3
		float vphase=phase-quadrant*(PI/2);
		sin_index=(vphase/(PI/2))*table_data.table_size; 
		cos_index=table_data.table_size-1-sin_index;
		if(quadrant&1) //in quadrant 1 and 3
		{
			temp_index=sin_index;
			sin_index=cos_index;
			cos_index=temp_index;
		}
		sin_sign=(quadrant>1)?-1:1; //in quadrant 2 and 3
		cos_sign=(quadrant&&quadrant<3)?-1:1; //in quadrant 1 and 2
		sinval=sin_sign*table_data.table[sin_index];
		cosval=cos_sign*table_data.table[cos_index];
		//we multiply two complex numbers.
		//how? enter this to maxima (software) for explanation:
		//   (a+b*%i)*(c+d*%i), rectform;
		iof(output,i)=cosval*iof(input,i)-sinval*qof(input,i);
		qof(output,i)=sinval*iof(input,i)+cosval*qof(input,i);
		phase+=phase_increment;
		while(phase>2*PI) phase-=2*PI; //@shift_math_cc: normalize phase
		while(phase<0) phase+=2*PI;
	}
	return phase;
}

int fir_decimate_cc(complexf *input, complexf *output, int input_size, int decimation, float *taps, int taps_length)
{
	//Theory: http://www.dspguru.com/dsp/faqs/multirate/decimation
	//It uses real taps. It returns the number of output samples actually written.
	//It needs overlapping input based on its returned value:
	//number of processed input samples = returned value * decimation factor
	//The output buffer should be at least input_length / 3.
	// i: input index | ti: tap index | oi: output index
	int oi=0;
	for(int i=0; i<input_size; i+=decimation) //@fir_decimate_cc: outer loop
	{
		if(i+taps_length>input_size) break;
		float acci=0;
		for(int ti=0; ti<taps_length; ti++) acci += (iof(input,i+ti)) * taps[ti]; //@fir_decimate_cc: i loop
		float accq=0;
		for(int ti=0; ti<taps_length; ti++) accq += (qof(input,i+ti)) * taps[ti]; //@fir_decimate_cc: q loop
		iof(output,oi)=acci;
		qof(output,oi)=accq;
		oi++;
	}
	return oi;
}

rational_resampler_ff_t rational_resampler_ff(float *input, float *output, int input_size, int interpolation, int decimation, float *taps, int taps_length, int last_taps_delay)
{
	
	//Theory: http://www.dspguru.com/dsp/faqs/multirate/resampling
	//oi: output index, i: tap index
	int output_size=input_size*interpolation/decimation;
	int oi;
	int startingi, delayi;
	//fprintf(stderr,"rational_resampler_ff | interpolation = %d | decimation = %d\ntaps_length = %d | input_size = %d | output_size = %d | last_taps_delay = %d\n",interpolation,decimation,taps_length,input_size,output_size,last_taps_delay);
	for (oi=0; oi<output_size; oi++) //@rational_resampler_ff (outer loop)
	{
		float acc=0;
		startingi=(oi*decimation+interpolation-1-last_taps_delay)/interpolation; //index of first input item to apply FIR on
		//delayi=startingi*interpolation-oi*decimation; //delay on FIR taps
		delayi=(last_taps_delay+startingi*interpolation-oi*decimation)%interpolation; //delay on FIR taps
		if(startingi+taps_length/interpolation+1>input_size) break; //we can't compute the FIR filter to some input samples at the end
		//fprintf(stderr,"outer loop | oi = %d | startingi = %d | taps delay = %d\n",oi,startingi,delayi);
		for(int i=0; i<(taps_length-delayi)/interpolation; i++)	//@rational_resampler_ff (inner loop)
		{
			//fprintf(stderr,"inner loop | input index = %d | tap index = %d | acc = %g\n",startingi+ii,i,acc);
			acc+=input[startingi+i]*taps[delayi+i*interpolation];
		}
		output[oi]=acc*interpolation;
	}
	rational_resampler_ff_t d;
	d.input_processed=startingi;
	d.output_size=oi;
	d.last_taps_delay=delayi;
	return d;
}

/*

The greatest challenge in resampling is figuring out which tap should be applied to which sample. 

Typical test patterns for rational_resampler_ff:

interpolation = 3, decimation = 1
values of [oi, startingi, taps delay] in the outer loop should be:
0 0 0
1 1 2
2 1 1
3 1 0
4 2 2
5 2 1

interpolation = 3, decimation = 2
values of [oi, startingi, taps delay] in the outer loop should be:
0 0 0
1 1 1
2 2 2
3 2 0
4 3 1
5 4 2

*/


void rational_resampler_get_lowpass_f(float* output, int output_size, int interpolation, int decimation, window_t window)
{

	//See 4.1.6 at: http://www.dspguru.com/dsp/faqs/multirate/resampling
	float cutoff_for_interpolation=1.0/interpolation;
	float cutoff_for_decimation=1.0/decimation;
	float cutoff = (cutoff_for_interpolation<cutoff_for_decimation)?cutoff_for_interpolation:cutoff_for_decimation; //get the lower
	firdes_lowpass_f(output, output_size, cutoff/2, window);
}

float inline fir_one_pass_ff(float* input, float* taps, int taps_length)
{
	float acc=0;
	for(int i=0;i<taps_length;i++) acc+=taps[i]*input[i]; //@fir_one_pass_ff
	return acc;
}

fractional_decimator_ff_t fractional_decimator_ff(float* input, float* output, int input_size, float rate, float *taps, int taps_length, fractional_decimator_ff_t d)
{
	if(rate<=1.0) return d; //sanity check, can't decimate <=1.0
	//This routine can handle floating point decimation rates. 
	//It linearly interpolates between two samples that are taken into consideration from the filtered input. 
	int oi=0;
	int index_high;
	float where=d.remain;
	float result_high, result_low;
	if(where==0.0) //in the first iteration index_high may be zero (so using the item index_high-1 would lead to invalid memory access). 
	{
		output[oi++]=fir_one_pass_ff(input,taps,taps_length);
		where+=rate;
	}

	int previous_index_high=-1;
	//we optimize to calculate ceilf(where) only once every iteration, so we do it here:
	for(;(index_high=ceilf(where))+taps_length<input_size;where+=rate) //@fractional_decimator_ff
	{
		if(previous_index_high==index_high-1) result_low=result_high; //if we step less than 2.0 then we do already have the result for the FIR filter for that index
		else result_low=fir_one_pass_ff(input+index_high-1,taps,taps_length); 
		result_high=fir_one_pass_ff(input+index_high,taps,taps_length);
		float register rate_between_samples=where-index_high+1;
		output[oi++]=result_low*(1-rate_between_samples)+result_high*rate_between_samples;
		previous_index_high=index_high;
	}

	d.input_processed=index_high-1;
	d.remain=where-d.input_processed;
	d.output_size=oi;
	return d;
}


void apply_fir_fft_cc(FFT_PLAN_T* plan, FFT_PLAN_T* plan_inverse, complexf* taps_fft, complexf* last_overlap, int overlap_size)
{
	//use the overlap & add method for filtering

	//calculate FFT on input buffer
	fft_execute(plan);

	//multiply the filter and the input
	complexf* in = plan->output;
	complexf* out = plan_inverse->input;
	
	for(int i=0;i<plan->size;i++) //@apply_fir_fft_cc: multiplication
	{
		iof(out,i)=iof(in,i)*iof(taps_fft,i)-qof(in,i)*qof(taps_fft,i);
		qof(out,i)=iof(in,i)*qof(taps_fft,i)+qof(in,i)*iof(taps_fft,i);
	}
	
	//calculate inverse FFT on multiplied buffer
	fft_execute(plan_inverse);
	
	//add the overlap of the previous segment
	complexf* result = plan_inverse->output;

	for(int i=0;i<plan->size;i++) //@apply_fir_fft_cc: normalize by fft_size
	{
		iof(result,i)/=plan->size;
		qof(result,i)/=plan->size;
	}
	
	for(int i=0;i<overlap_size;i++) //@apply_fir_fft_cc: add overlap
	{
		iof(result,i)=iof(result,i)+iof(last_overlap,i);
		qof(result,i)=qof(result,i)+qof(last_overlap,i);
	}
	
}

/*
           __  __       _                          _       _       _                 
     /\   |  \/  |     | |                        | |     | |     | |                
    /  \  | \  / |   __| | ___ _ __ ___   ___   __| |_   _| | __ _| |_ ___  _ __ ___ 
   / /\ \ | |\/| |  / _` |/ _ \ '_ ` _ \ / _ \ / _` | | | | |/ _` | __/ _ \| '__/ __|
  / ____ \| |  | | | (_| |  __/ | | | | | (_) | (_| | |_| | | (_| | || (_) | |  \__ \
 /_/    \_\_|  |_|  \__,_|\___|_| |_| |_|\___/ \__,_|\__,_|_|\__,_|\__\___/|_|  |___/
                                                                                                                                                                        
*/

void amdemod_cf(complexf* input, float *output, int input_size)
{
	//@amdemod: i*i+q*q
	for (int i=0; i<input_size; i++)
	{ 
		output[i]=iof(input,i)*iof(input,i)+qof(input,i)*qof(input,i);
	}
	//@amdemod: sqrt
	for (int i=0; i<input_size; i++)
	{ 
		output[i]=sqrt(output[i]);
	}
}

void amdemod_estimator_cf(complexf* input, float *output, int input_size, float alpha, float beta)
{
	//concept is explained here:
	//http://www.dspguru.com/dsp/tricks/magnitude-estimator

	//default: optimize for min RMS error
	if(alpha==0) 
	{
		alpha=0.947543636291;
		beta=0.392485425092;
	}

	//@amdemod_estimator
	for (int i=0; i<input_size; i++)
	{ 
		float abs_i=iof(input,i);
		if(abs_i<0) abs_i=-abs_i;
		float abs_q=qof(input,i);
		if(abs_q<0) abs_q=-abs_q;
		float max_iq=abs_i;
		if(abs_q>max_iq) max_iq=abs_q;
		float min_iq=abs_i;
		if(abs_q<min_iq) min_iq=abs_q;

		output[i]=alpha*max_iq+beta*min_iq;
	}
}

dcblock_preserve_t dcblock_ff(float* input, float* output, int input_size, float a, dcblock_preserve_t preserved)
{
	//after AM demodulation, a DC blocking filter should be used to remove the DC component from the signal.
	//Concept: http://peabody.sapp.org/class/dmp2/lab/dcblock/
	//output size equals to input_size;
	//preserve can be initialized to zero on first run.
	if(a==0) a=0.999; //default value, simulate in octave: freqz([1 -1],[1 -0.99])
	output[0]=input[0]-preserved.last_input+a*preserved.last_output;
	for(int i=1; i<input_size; i++) //@dcblock_f
	{
		output[i]=input[i]-input[i-1]+a*output[i-1];
	}
	preserved.last_input=input[input_size-1];
	preserved.last_output=output[input_size-1];
	return preserved;
}

float fastdcblock_ff(float* input, float* output, int input_size, float last_dc_level)
{
	//this DC block filter does moving average block-by-block.
	//this is the most computationally efficient
	//input and output buffer is allowed to be the same
	//http://www.digitalsignallabs.com/dcblock.pdf
	float avg=0.0;
	for(int i=0;i<input_size;i++) //@fastdcblock_ff: calculate block average
	{
		avg+=input[i];
	}
	avg/=input_size;
	
	float avgdiff=avg-last_dc_level; 
	//DC removal level will change lineraly from last_dc_level to avg.
	for(int i=0;i<input_size;i++) //@fastdcblock_ff: remove DC component
	{
		float dc_removal_level=last_dc_level+avgdiff*((float)i/input_size);
		output[i]=input[i]-dc_removal_level;
	}
	return avg;
}

#define FASTAGC_MAX_GAIN (65e3)

void fastagc_ff(fastagc_ff_t* input, float* output)
{
	//Gain is processed on blocks of samples.
	//You have to supply three blocks of samples before the first block comes out.
	//AGC reaction speed equals input_size*samp_rate*2

	//The algorithm calculates target gain at the end of the first block out of the peak value of all the three blocks. 
	//This way the gain change can easily react if there is any peak in the third block.
	//Pros: can be easily speeded up with loop vectorization, easy to implement
	//Cons: needs 3 buffers, dos not behave similarly to real AGC circuits

	//Get the peak value of new input buffer
	float peak_input=0;
	for(int i=0;i<input->input_size;i++) //@fastagc_ff: peak search
	{
		float val=fabs(input->buffer_input[i]);
		if(val>peak_input) peak_input=val;
	}

	//Determine the maximal peak out of the three blocks
	float target_peak=peak_input;
	if(target_peak<input->peak_2) target_peak=input->peak_2;
	if(target_peak<input->peak_1) target_peak=input->peak_1;
	
	//we change the gain linearly on the apply_block from the last_gain to target_gain.
	float target_gain=input->reference/target_peak;
	if(target_gain>FASTAGC_MAX_GAIN) target_gain=FASTAGC_MAX_GAIN;

	for(int i=0;i<input->input_size;i++) //@fastagc_ff: apply gain
	{
		float rate=(float)i/input->input_size;
		float gain=input->last_gain*(1.0-rate)+target_gain*rate;
		output[i]=input->buffer_1[i]*gain;
	}

	//Shift the three buffers
	float* temp_pointer=input->buffer_1;
	input->buffer_1=input->buffer_2;
	input->peak_1=input->peak_2;
	input->buffer_2=input->buffer_input;
	input->peak_2=peak_input;
	input->buffer_input=temp_pointer;
	input->last_gain=target_gain;
	//fprintf(stderr,"target_gain=%g\n", target_gain);
}


/*
  ______ __  __        _                          _       _       _                 
 |  ____|  \/  |      | |                        | |     | |     | |                
 | |__  | \  / |    __| | ___ _ __ ___   ___   __| |_   _| | __ _| |_ ___  _ __ ___ 
 |  __| | |\/| |   / _` |/ _ \ '_ ` _ \ / _ \ / _` | | | | |/ _` | __/ _ \| '__/ __|
 | |    | |  | |  | (_| |  __/ | | | | | (_) | (_| | |_| | | (_| | || (_) | |  \__ \
 |_|    |_|  |_|   \__,_|\___|_| |_| |_|\___/ \__,_|\__,_|_|\__,_|\__\___/|_|  |___/

*/


float fmdemod_atan_cf(complexf* input, float *output, int input_size, float last_phase)
{
	//GCC most likely won't vectorize nor atan, nor atan2.
	//For more comments, look at: https://github.com/simonyiszk/minidemod/blob/master/minidemod-wfm-atan.c
	float phase, dphase;
	for (int i=0; i<input_size; i++) //@fmdemod_atan_novect
	{
		phase=argof(input,i);
		dphase=phase-last_phase;
		if(dphase<-PI) dphase+=2*PI;
		if(dphase>PI) dphase-=2*PI;
		output[i]=dphase/PI;
		last_phase=phase;
	}
	return last_phase;
}

#define fmdemod_quadri_K 0.340447550238101026565118445432744920253753662109375
//this constant ensures proper scaling for qa_fmemod testcases for SNR calculation and more.

complexf fmdemod_quadri_novect_cf(complexf* input, float* output, int input_size, complexf last_sample)
{
	output[0]=fmdemod_quadri_K*(iof(input,0)*(qof(input,0)-last_sample.q)-qof(input,0)*(iof(input,0)-last_sample.i))/(iof(input,0)*iof(input,0)+qof(input,0)*qof(input,0));
	for (int i=1; i<input_size; i++) //@fmdemod_quadri_novect_cf
	{
		float qnow=qof(input,i);
		float qlast=qof(input,i-1);
		float inow=iof(input,i);
		float ilast=iof(input,i-1);
		output[i]=fmdemod_quadri_K*(inow*(qnow-qlast)-qnow*(inow-ilast))/(inow*inow+qnow*qnow);
		//TODO: expression can be simplified as: (qnow*ilast-inow*qlast)/(inow*inow+qnow*qnow)
	}
	return input[input_size-1];
}


complexf fmdemod_quadri_cf(complexf* input, float* output, int input_size, float *temp, complexf last_sample)
{
	float* temp_dq=temp;
	float* temp_di=temp+input_size;

	temp_dq[0]=qof(input,0)-last_sample.q;
	for (int i=1; i<input_size; i++) //@fmdemod_quadri_cf: dq
	{		
		temp_dq[i]=qof(input,i)-qof(input,i-1);
	}

	temp_di[0]=iof(input,0)-last_sample.i;
	for (int i=1; i<input_size; i++) //@fmdemod_quadri_cf: di
	{		
		temp_di[i]=iof(input,i)-iof(input,i-1);
	}

	for (int i=0; i<input_size; i++) //@fmdemod_quadri_cf: output numerator
	{		
		output[i]=(iof(input,i)*temp_dq[i]-qof(input,i)*temp_di[i]);
	}
	for (int i=0; i<input_size; i++) //@fmdemod_quadri_cf: output denomiator
	{		
		temp[i]=iof(input,i)*iof(input,i)+qof(input,i)*qof(input,i);
	}	
	for (int i=0; i<input_size; i++) //@fmdemod_quadri_cf: output division
	{		
		output[i]=fmdemod_quadri_K*output[i]/temp[i];
	}
	
	return input[input_size-1];
}

inline int is_nan(float f) 
{
	//http://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c
    unsigned u = *(unsigned*)&f;
    return (u&0x7F800000) == 0x7F800000 && (u&0x7FFFFF); // Both NaN and qNan.
}


float deemphasis_wfm_ff (float* input, float* output, int input_size, float tau, int sample_rate, float last_output)
{
	/* 
		typical time constant (tau) values:
		WFM transmission in USA: 75 us -> tau = 75e-6
		WFM transmission in EU:  50 us -> tau = 50e-6
		More info at: http://www.cliftonlaboratories.com/fm_receivers_and_de-emphasis.htm
		Simulate in octave: tau=75e-6; dt=1/48000; alpha = dt/(tau+dt); freqz([alpha],[1 -(1-alpha)])
	*/
	float dt = 1.0/sample_rate;
	float alpha = dt/(tau+dt);
	if(is_nan(last_output)) last_output=0.0; //if last_output is NaN
	output[0]=alpha*input[0]+(1-alpha)*last_output; 
	for (int i=1;i<input_size;i++) //@deemphasis_wfm_ff
       output[i]=alpha*input[i]+(1-alpha)*output[i-1]; //this is the simplest IIR LPF
   	return output[input_size-1];
}

#define DNFMFF_ADD_ARRAY(x) if(sample_rate==x) { taps=deemphasis_nfm_predefined_fir_##x; taps_length=sizeof(deemphasis_nfm_predefined_fir_##x)/sizeof(float); }

int deemphasis_nfm_ff (float* input, float* output, int input_size, int sample_rate)
{
	/*
		Warning! This only works on predefined samplerates, as it uses fixed FIR coefficients defined in predefined.h
		However, there are the octave commands to generate the taps for your custom (fixed) sample rate.
		What it does: 
			- reject below 400 Hz
			- passband between 400 HZ - 4 kHz, but with 20 dB/decade rolloff (for deemphasis)
			- reject everything above 4 kHz
	*/
	float* taps;
	int taps_length=0;

	DNFMFF_ADD_ARRAY(48000)
	DNFMFF_ADD_ARRAY(44100)
	DNFMFF_ADD_ARRAY(8000)
	DNFMFF_ADD_ARRAY(11025)

	if(!taps_length) return 0; //sample rate n
	int i;
	for(i=0;i<input_size-taps_length;i++) //@deemphasis_nfm_ff: outer loop
	{
		float acc=0;
		for(int ti=0;ti<taps_length;ti++) acc+=taps[ti]*input[i+ti]; //@deemphasis_nfm_ff: inner loop
		output[i]=acc;
	}
	return i; //number of samples processed (and output samples)
}

void limit_ff(float* input, float* output, int input_size, float max_amplitude)
{
	for (int i=0; i<input_size; i++) //@limit_ff
	{		
		output[i]=(max_amplitude<input[i])?max_amplitude:input[i];
		output[i]=(-max_amplitude>output[i])?-max_amplitude:output[i];
	}	
}

void gain_ff(float* input, float* output, int input_size, float gain)
{
	for(int i=0;i<input_size;i++) output[i]=gain*input[i]; //@gain_ff
}


/*
  ______        _     ______               _             _______                   __                     
 |  ____|      | |   |  ____|             (_)           |__   __|                 / _|                    
 | |__ __ _ ___| |_  | |__ ___  _   _ _ __ _  ___ _ __     | |_ __ __ _ _ __  ___| |_ ___  _ __ _ __ ___  
 |  __/ _` / __| __| |  __/ _ \| | | | '__| |/ _ \ '__|    | | '__/ _` | '_ \/ __|  _/ _ \| '__| '_ ` _ \ 
 | | | (_| \__ \ |_  | | | (_) | |_| | |  | |  __/ |       | | | | (_| | | | \__ \ || (_) | |  | | | | | |
 |_|  \__,_|___/\__| |_|  \___/ \__,_|_|  |_|\___|_|       |_|_|  \__,_|_| |_|___/_| \___/|_|  |_| |_| |_|                                                                                                                                                                                                                   

*/

int log2n(int x)
{
	int result=-1;
	for(int i=0;i<31;i++) 
	{
		if((x>>i)&1) //@@log2n
		{
			if (result==-1) result=i;
			else return -1;
		}
	}
	return result;
}

int next_pow2(int x)
{
	int pow2;
	//portability? (31 is the problem)
	for(int i=0;i<31;i++) 
	{
		if(x<(pow2=1<<i)) return pow2; //@@next_pow2
	}
	return -1;
}

void apply_window_c(complexf* input, complexf* output, int size, window_t window)
{
	float (*window_function)(float)=firdes_get_window_kernel(window);
	for(int i=0;i<size;i++) //@apply_window_c
	{
		float rate=(float)i/(size-1);
		iof(output,i)=iof(input,i)*window_function(2.0*rate+1.0);
		qof(output,i)=qof(input,i)*window_function(2.0*rate+1.0);
	}
}

void apply_window_f(float* input, float* output, int size, window_t window)
{
	float (*window_function)(float)=firdes_get_window_kernel(window);
	for(int i=0;i<size;i++) //@apply_window_f
	{
		float rate=(float)i/(size-1);
		output[i]=input[i]*window_function(2.0*rate+1.0);
	}
}

void logpower_cf(complexf* input, float* output, int size, float add_db)
{
	for(int i=0;i<size;i++) output[i]=iof(input,i)*iof(input,i) + qof(input,i)*qof(input,i); //@logpower_cf: pass 1
	
	for(int i=0;i<size;i++) output[i]=log10(output[i]); //@logpower_cf: pass 2

	for(int i=0;i<size;i++) output[i]=10*output[i]+add_db; //@logpower_cf: pass 3
}

/*
  _____        _                                            _             
 |  __ \      | |                                          (_)            
 | |  | | __ _| |_ __ _    ___ ___  _ ____   _____ _ __ ___ _  ___  _ __  
 | |  | |/ _` | __/ _` |  / __/ _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \ 
 | |__| | (_| | || (_| | | (_| (_) | | | \ V /  __/ |  \__ \ | (_) | | | |
 |_____/ \__,_|\__\__,_|  \___\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
                                                                          
*/                                                                          

void convert_u8_f(unsigned char* input, float* output, int input_size)
{
	for(int i=0;i<input_size;i++) output[i]=((float)input[i])/(UCHAR_MAX/2.0)-1.0; //@convert_u8_f
}

void convert_i8_f(signed char* input, float* output, int input_size)
{
	for(int i=0;i<input_size;i++) output[i]=((float)input[i])/SCHAR_MAX; //@convert_i8_f
}

void convert_i16_f(short* input, float* output, int input_size)
{
	for(int i=0;i<input_size;i++) output[i]=(float)input[i]/SHRT_MAX; //@convert_i16_f
}

void convert_f_u8(float* input, unsigned char* output, int input_size)
{
	for(int i=0;i<input_size;i++) output[i]=input[i]*UCHAR_MAX*0.5+128; //@convert_f_u8
	//128 above is the correct value to add. In any other case a DC component 
	//of at least -60 dB is shown on the FFT plot after convert_f_u8 -> convert_u8_f
}

void convert_f_i8(float* input, signed char* output, int input_size)
{
	/*for(int i=0;i<input_size;i++)
	{
		if(input[i]>1.0) input[i]=1.0;
		if(input[i]<-1.0) input[i]=-1.0;
	}*/
	for(int i=0;i<input_size;i++) output[i]=input[i]*SCHAR_MAX; //@convert_f_i8
}
void convert_f_i16(float* input, short* output, int input_size)
{
	/*for(int i=0;i<input_size;i++) 
	{
		if(input[i]>1.0) input[i]=1.0;
		if(input[i]<-1.0) input[i]=-1.0;
	}*/
	for(int i=0;i<input_size;i++) output[i]=input[i]*SHRT_MAX; //@convert_f_i16
}

int trivial_vectorize()
{
	//this function is trivial to vectorize and should pass on both NEON and SSE
	int a[1024], b[1024], c[1024];
	for(int i=0; i<1024; i++) //@trivial_vectorize: should pass :-)
	{
		c[i]=a[i]*b[i];
	}
	return c[0];
}


