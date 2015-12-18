/*
This software is part of libcsdr, a set of simple DSP routines for 
Software Defined Radio.

Copyright (c) 2014-2015, Andras Retzler <randras@sdr.hu>
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
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>

#include "libcsdr.h"
#include "libcsdr_gpl.h"

#define T_BUFSIZE (1024*1024/4)
#define T_N (200)
#define T_TAPS (1023)
#define T_DECFACT (200)

int main()
{
	fprintf(stderr,"Getting a %d of random samples...\n", T_BUFSIZE);
	int urand_fp = open("/dev/urandom",O_RDWR);
	unsigned char* buf_u8 = (unsigned char*)malloc(sizeof(unsigned char)*T_BUFSIZE*2);
	complexf* buf_c = (complexf*)malloc(sizeof(complexf)*T_BUFSIZE);
	complexf* outbuf_c = (complexf*)malloc(sizeof(complexf)*T_BUFSIZE);
	read(urand_fp, buf_u8, T_BUFSIZE);
	close(urand_fp);
	
	for(int i=0;i<T_BUFSIZE;i++)
	{ 
		iof(buf_c,i)=buf_u8[2*i]/128.0;
		qof(buf_c,i)=buf_u8[2*i+1]/128.0;
	}


	float* taps_f = (float*)malloc(sizeof(float)*T_TAPS);
	firdes_lowpass_f(taps_f, T_TAPS, 1.0f/T_DECFACT, WINDOW_DEFAULT);

	struct timespec start_time, end_time;	

	fprintf(stderr,"Starting tests of processing %d samples...\n", T_BUFSIZE*T_N);

	//fir_decimate_cc
        clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
        for(int i=0;i<T_N;i++) fir_decimate_cc(buf_c, outbuf_c, T_BUFSIZE, 10, taps_f, T_TAPS);
        clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);
        fprintf(stderr,"fir_decimate_cc done in %g seconds.\n",TIME_TAKEN(start_time,end_time));


	//shift_math_cc
	float starting_phase = 0;

	clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
	for(int i=0;i<T_N;i++) starting_phase = shift_math_cc(buf_c, outbuf_c, T_BUFSIZE, 0.1, starting_phase);
	clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);
	fprintf(stderr,"shift_math_cc done in %g seconds.\n",TIME_TAKEN(start_time,end_time));

	//shift_table_cc	
	shift_table_data_t shift_table_data=shift_table_init(65536);
	starting_phase = 0;

	clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
	for(int i=0;i<T_N;i++) starting_phase = starting_phase=shift_table_cc(buf_c, outbuf_c, T_BUFSIZE, 0.1, shift_table_data, starting_phase);;
	clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);
	fprintf(stderr,"shift_table_cc (table size = %d) done in %g seconds.\n",65536,TIME_TAKEN(start_time,end_time));


	//shift_addition_cc	
	shift_addition_data_t data_addition = shift_addition_init(0.1);
	starting_phase = 0;

	clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
	for(int i=0;i<T_N;i++) starting_phase = shift_addition_cc(buf_c, outbuf_c, T_BUFSIZE, data_addition, starting_phase);
	clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);
	fprintf(stderr,"shift_addition_cc done in %g seconds.\n",TIME_TAKEN(start_time,end_time));

	//shift_addfast_cc	
	shift_addfast_data_t data_addfast = shift_addfast_init(0.1);
	starting_phase = 0;

	clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
	for(int i=0;i<T_N;i++) starting_phase = shift_addfast_cc(buf_c, outbuf_c, T_BUFSIZE, &data_addfast, starting_phase);
	clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);
	fprintf(stderr,"shift_addfast_cc done in %g seconds.\n",TIME_TAKEN(start_time,end_time));

	//shift_unroll_cc
	shift_unroll_data_t data_unroll = shift_unroll_init(0.1, T_BUFSIZE);
	starting_phase = 0;

	clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
	for(int i=0;i<T_N;i++) starting_phase = shift_unroll_cc(buf_c, outbuf_c, T_BUFSIZE, &data_unroll, starting_phase);
	clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);
	fprintf(stderr,"shift_unroll_cc done in %g seconds.\n",TIME_TAKEN(start_time,end_time));


}
