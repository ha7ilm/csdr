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

#include "fastddc.h"

//DDC implementation based on:
//http://www.3db-labs.com/01598092_MultibandFilterbank.pdf

inline int is_integer(float a) { return floorf(a) == a; }

int fastddc_init(fastddc_t* ddc, float transition_bw, int decimation, float shift_rate)
{
	ddc->pre_decimation = 1; //this will be done in the frequency domain
	ddc->post_decimation = decimation; //this will be done in the time domain
	while( is_integer((float)ddc->post_decimation/2) && ddc->post_decimation/2 != 1) 
	{
		ddc->post_decimation/=2;
		ddc->pre_decimation*=2;
	}
	ddc->taps_real_length = firdes_filter_len(transition_bw); //the number of non-zero taps
	ddc->taps_length = ceil(ddc->taps_real_length/(float)ddc->pre_decimation) * ddc->pre_decimation; //the number of taps must be a multiple of the decimation factor
	ddc->fft_size = next_pow2(ddc->taps_length * 4); //it is a good rule of thumb for performance (based on the article), but we should do benchmarks
	while (ddc->fft_size<ddc->pre_decimation) ddc->fft_size*=2; //fft_size should be a multiple of pre_decimation
	ddc->overlap_length = ddc->taps_length - 1;
	ddc->input_size = ddc->fft_size - ddc->overlap_length;
	ddc->fft_inv_size = ddc->fft_size / ddc->pre_decimation;

	//Shift operation in the frequency domain: we can shift by a multiple of v.
	ddc->v = ddc->fft_size/ddc->overlap_length; //+-1 ? (or maybe ceil() this?) //TODO: why?
	int middlebin=ddc->fft_size / 2;
	ddc->startbin = middlebin + middlebin * shift_rate * 2;	
	//fprintf(stderr, "ddc->startbin=%g\n",(float)ddc->startbin);
	ddc->startbin = ddc->v * round( ddc->startbin / (float)ddc->v );
	//fprintf(stderr, "ddc->startbin=%g\n",(float)ddc->startbin);
	ddc->offsetbin = ddc->startbin - middlebin;
	ddc->post_shift = shift_rate-((float)ddc->offsetbin/ddc->fft_size);
	ddc->pre_shift = ddc->offsetbin/(float)ddc->fft_size;
	ddc->dsadata = decimating_shift_addition_init(ddc->post_shift, ddc->post_decimation);

	//Overlap is scrapd, not added
	ddc->scrap=ddc->overlap_length/ddc->pre_decimation; //TODO this is problematic sometimes! overlap_length = 401 :: scrap = 200
	ddc->post_input_size=ddc->fft_inv_size-ddc->scrap;

	return ddc->fft_size<=2; //returns true on error
}


void fastddc_print(fastddc_t* ddc, char* source)
{
	fprintf(stderr,
		"%s: fastddc_print_sizes(): (fft_size = %d) = (taps_length = %d) + (input_size = %d) - 1\n"
		"  overlap     ::  (overlap_length = %d) = taps_length - 1, taps_real_length = %d\n"
		"  decimation  ::  decimation = (pre_decimation = %d) * (post_decimation = %d), fft_inv_size = %d\n"
		"  shift       ::  startbin = %d, offsetbin = %d, v = %d, pre_shift = %g, post_shift = %g\n"
		"  o&s         ::  post_input_size = %d, scrap = %d\n"
		, 
		source, ddc->fft_size, ddc->taps_length, ddc->input_size, 
		ddc->overlap_length, ddc->taps_real_length,
		ddc->pre_decimation, ddc->post_decimation, ddc->fft_inv_size,
		ddc->startbin, ddc->offsetbin, ddc->v, ddc->pre_shift, ddc->post_shift, 
		ddc->post_input_size, ddc->scrap );
}

void fft_swap_sides(complexf* io, int fft_size)
{
	int middle=fft_size/2;
	complexf temp;
	for(int i=0;i<middle;i++)
	{
		iof(&temp,0)=iof(io,i);
		qof(&temp,0)=qof(io,i);
		iof(io,i)=iof(io,i+middle);
		qof(io,i)=qof(io,i+middle);
		iof(io,i+middle)=iof(&temp,0);
		qof(io,i+middle)=qof(&temp,0);
	}
}

decimating_shift_addition_status_t fastddc_inv_cc(complexf* input, complexf* output, fastddc_t* ddc, FFT_PLAN_T* plan_inverse, complexf* taps_fft, decimating_shift_addition_status_t shift_stat)
{
	//implements DDC by using the overlap & scrap method
	//TODO: +/-1s on overlap_size et al
	//input shoud have ddc->fft_size number of elements

	complexf* inv_input = plan_inverse->input;
	complexf* inv_output = plan_inverse->output;

	//Initialize buffers for inverse FFT to zero
	for(int i=0;i<plan_inverse->size;i++)
	{
		iof(inv_input,i)=0;
		qof(inv_input,i)=0;
	}

	//Alias & shift & filter at once
	fft_swap_sides(input, ddc->fft_size); //this is not very optimal, but now we stick with this slow solution until we got the algorithm working
	for(int i=0;i<ddc->fft_size;i++)
	{
		int output_index = (ddc->startbin+i)%plan_inverse->size;
		int tap_index = (ddc->fft_size+i-ddc->offsetbin)%ddc->fft_size;
		cmultadd(inv_input+output_index, input+i, taps_fft+tap_index); //cmultadd(output, input1, input2):   complex output += complex input1 * complex input 2
	}

	fft_execute(plan_inverse);
	fft_swap_sides(inv_output,plan_inverse->size);

	//Normalize data
	for(int i=0;i<plan_inverse->size;i++) //@fastddc_inv_cc: normalize by size
	{
		iof(inv_output,i)/=plan_inverse->size;
		qof(inv_output,i)/=plan_inverse->size;
	}
	
	//Overlap is scrapped, not added
	//Shift correction
	shift_stat=decimating_shift_addition_cc(inv_output+ddc->scrap, output, ddc->post_input_size, ddc->dsadata, ddc->post_decimation, shift_stat);
	//memcpy(inv_output+ddc->scrap, output
	return shift_stat;
}
