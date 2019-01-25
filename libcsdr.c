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
#include <stdarg.h>

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
{   //"Dummy" window kernel, do not use; an unwindowed FIR filter may have bad frequency response
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

void normalize_fir_f(float* input, float* output, int length)
{
    //Normalize filter kernel
    float sum=0;
    for(int i=0;i<length;i++) //@normalize_fir_f: normalize pass 1
        sum+=input[i];
    for(int i=0;i<length;i++) //@normalize_fir_f: normalize pass 2
        output[i]=input[i]/sum;
}

void firdes_lowpass_f(float *output, int length, float cutoff_rate, window_t window)
{   //Generates symmetric windowed sinc FIR filter real taps
    //  length should be odd
    //  cutoff_rate is (cutoff frequency/sampling frequency)
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
    normalize_fir_f(output,output,length);
}

void firdes_bandpass_c(complexf *output, int length, float lowcut, float highcut, window_t window)
{
    //To generate a complex filter:
    //  1. we generate a real lowpass filter with a bandwidth of highcut-lowcut
    //  2. we shift the filter taps spectrally by multiplying with e^(j*w), so we get complex taps
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


shift_unroll_data_t shift_unroll_init(float rate, int size)
{
    shift_unroll_data_t output;
    output.phase_increment=2*rate*PI;
    output.size = size;
    output.dsin=(float*)malloc(sizeof(float)*size);
    output.dcos=(float*)malloc(sizeof(float)*size);
    float myphase = 0;
    for(int i=0;i<size;i++)
    {
        myphase += output.phase_increment;
        while(myphase>PI) myphase-=2*PI;
        while(myphase<-PI) myphase+=2*PI;
        output.dsin[i]=sin(myphase);
        output.dcos[i]=cos(myphase);
    }
    return output;
}

float shift_unroll_cc(complexf *input, complexf* output, int input_size, shift_unroll_data_t* d, float starting_phase)
{
    //input_size should be multiple of 4
    //fprintf(stderr, "shift_addfast_cc: input_size = %d\n", input_size);
    float cos_start=cos(starting_phase);
    float sin_start=sin(starting_phase);
    register float cos_val, sin_val;
    for(int i=0;i<input_size; i++) //@shift_unroll_cc
    {
        cos_val = cos_start * d->dcos[i] - sin_start * d->dsin[i];
        sin_val  = sin_start * d->dcos[i] + cos_start * d->dsin[i];
        iof(output,i)=cos_val*iof(input,i)-sin_val*qof(input,i);
        qof(output,i)=sin_val*iof(input,i)+cos_val*qof(input,i);
    }
    starting_phase+=input_size*d->phase_increment;
    while(starting_phase>PI) starting_phase-=2*PI;
    while(starting_phase<-PI) starting_phase+=2*PI;
    return starting_phase;
}

shift_addfast_data_t shift_addfast_init(float rate)
{
    shift_addfast_data_t output;
    output.phase_increment=2*rate*PI;
    for(int i=0;i<4;i++)
    {
        output.dsin[i]=sin(output.phase_increment*(i+1));
        output.dcos[i]=cos(output.phase_increment*(i+1));
    }
    return output;
}

#ifdef NEON_OPTS
#pragma message "Manual NEON optimizations are ON: we have a faster shift_addfast_cc now."

float shift_addfast_cc(complexf *input, complexf* output, int input_size, shift_addfast_data_t* d, float starting_phase)
{
    //input_size should be multiple of 4
    float cos_start[4], sin_start[4];
    float cos_vals[4], sin_vals[4];
    for(int i=0;i<4;i++)
    {
        cos_start[i] = cos(starting_phase);
        sin_start[i] = sin(starting_phase);
    }

    float* pdcos = d->dcos;
    float* pdsin = d->dsin;
    register float* pinput = (float*)input;
    register float* pinput_end = (float*)(input+input_size);
    register float* poutput = (float*)output;

    //Register map:
    #define RDCOS "q0" //dcos, dsin
    #define RDSIN "q1"
    #define RCOSST "q2" //cos_start, sin_start
    #define RSINST "q3"
    #define RCOSV "q4" //cos_vals, sin_vals
    #define RSINV "q5"
    #define ROUTI "q6" //output_i, output_q
    #define ROUTQ "q7"
    #define RINPI "q8" //input_i, input_q
    #define RINPQ "q9"
    #define R3(x,y,z) x ", " y ", " z "\n\t"

    asm volatile( //(the range of q is q0-q15)
        "       vld1.32 {" RDCOS "}, [%[pdcos]]\n\t"
        "       vld1.32 {" RDSIN "}, [%[pdsin]]\n\t"
        "       vld1.32 {" RCOSST "}, [%[cos_start]]\n\t"
        "       vld1.32 {" RSINST "}, [%[sin_start]]\n\t"
        "for_addfast: vld2.32 {" RINPI "-" RINPQ "}, [%[pinput]]!\n\t" //load q0 and q1 directly from the memory address stored in pinput, with interleaving (so that we get the I samples in RINPI and the Q samples in RINPQ), also increment the memory address in pinput (hence the "!" mark)

        //C version:
        //cos_vals[j] = cos_start * d->dcos[j] - sin_start * d->dsin[j];
        //sin_vals[j] = sin_start * d->dcos[j] + cos_start * d->dsin[j];

        "       vmul.f32 " R3(RCOSV, RCOSST, RDCOS)  //cos_vals[i] = cos_start * d->dcos[i]
        "       vmls.f32 " R3(RCOSV, RSINST, RDSIN)  //cos_vals[i] -= sin_start * d->dsin[i]
        "       vmul.f32 " R3(RSINV, RSINST, RDCOS)  //sin_vals[i] = sin_start * d->dcos[i]
        "       vmla.f32 " R3(RSINV, RCOSST, RDSIN)  //sin_vals[i] += cos_start * d->dsin[i]

        //C version:
        //iof(output,4*i+j)=cos_vals[j]*iof(input,4*i+j)-sin_vals[j]*qof(input,4*i+j);
        //qof(output,4*i+j)=sin_vals[j]*iof(input,4*i+j)+cos_vals[j]*qof(input,4*i+j);
        "       vmul.f32 " R3(ROUTI, RCOSV, RINPI) //output_i =  cos_vals * input_i
        "       vmls.f32 " R3(ROUTI, RSINV, RINPQ) //output_i -= sin_vals * input_q
        "       vmul.f32 " R3(ROUTQ, RSINV, RINPI) //output_q =  sin_vals * input_i
        "       vmla.f32 " R3(ROUTQ, RCOSV, RINPQ) //output_i += cos_vals * input_q

        "       vst2.32 {" ROUTI "-" ROUTQ "}, [%[poutput]]!\n\t" //store the outputs in memory
        //"     add %[poutput],%[poutput],#32\n\t"
        "       vdup.32 " RCOSST ", d9[1]\n\t" // cos_start[0-3] = cos_vals[3]
        "       vdup.32 " RSINST ", d11[1]\n\t" // sin_start[0-3] = sin_vals[3]

        "       cmp %[pinput], %[pinput_end]\n\t" //if(pinput != pinput_end)
        "       bcc for_addfast\n\t"              //    then goto for_addfast
    :
        [pinput]"+r"(pinput), [poutput]"+r"(poutput) //output operand list -> C variables that we will change from ASM
    :
        [pinput_end]"r"(pinput_end), [pdcos]"r"(pdcos), [pdsin]"r"(pdsin), [sin_start]"r"(sin_start), [cos_start]"r"(cos_start) //input operand list
    :
        "memory", "q0", "q1", "q2", "q4", "q5", "q6", "q7", "q8", "q9", "cc" //clobber list
    );
    starting_phase+=input_size*d->phase_increment;
    while(starting_phase>PI) starting_phase-=2*PI;
    while(starting_phase<-PI) starting_phase+=2*PI;
    return starting_phase;
}

#else


#if 1

#define SADF_L1(j) cos_vals_ ## j = cos_start * dcos_ ## j - sin_start * dsin_ ## j; \
    sin_vals_ ## j = sin_start * dcos_ ## j + cos_start * dsin_ ## j;
#define SADF_L2(j) iof(output,4*i+j)=(cos_vals_ ## j)*iof(input,4*i+j)-(sin_vals_ ## j)*qof(input,4*i+j); \
    qof(output,4*i+j)=(sin_vals_ ## j)*iof(input,4*i+j)+(cos_vals_ ## j)*qof(input,4*i+j);

float shift_addfast_cc(complexf *input, complexf* output, int input_size, shift_addfast_data_t* d, float starting_phase)
{
    //input_size should be multiple of 4
    //fprintf(stderr, "shift_addfast_cc: input_size = %d\n", input_size);
    float cos_start=cos(starting_phase);
    float sin_start=sin(starting_phase);
    float register cos_vals_0, cos_vals_1, cos_vals_2, cos_vals_3,
        sin_vals_0, sin_vals_1, sin_vals_2, sin_vals_3,
        dsin_0 = d->dsin[0], dsin_1 = d->dsin[1], dsin_2 = d->dsin[2], dsin_3 = d->dsin[3],
        dcos_0 = d->dcos[0], dcos_1 = d->dcos[1], dcos_2 = d->dcos[2], dcos_3 = d->dcos[3];

    for(int i=0;i<input_size/4; i++) //@shift_addfast_cc
    {
        SADF_L1(0)
        SADF_L1(1)
        SADF_L1(2)
        SADF_L1(3)
        SADF_L2(0)
        SADF_L2(1)
        SADF_L2(2)
        SADF_L2(3)
        cos_start = cos_vals_3;
        sin_start = sin_vals_3;
    }
    starting_phase+=input_size*d->phase_increment;
    while(starting_phase>PI) starting_phase-=2*PI;
    while(starting_phase<-PI) starting_phase+=2*PI;
    return starting_phase;
}
#else
float shift_addfast_cc(complexf *input, complexf* output, int input_size, shift_addfast_data_t* d, float starting_phase)
{
    //input_size should be multiple of 4
    //fprintf(stderr, "shift_addfast_cc: input_size = %d\n", input_size);
    float cos_start=cos(starting_phase);
    float sin_start=sin(starting_phase);
    float cos_vals[4], sin_vals[4];
    for(int i=0;i<input_size/4; i++) //@shift_addfast_cc
    {
        for(int j=0;j<4;j++) //@shift_addfast_cc
        {
            cos_vals[j] = cos_start * d->dcos[j] - sin_start * d->dsin[j];
            sin_vals[j] = sin_start * d->dcos[j] + cos_start * d->dsin[j];
        }
        for(int j=0;j<4;j++) //@shift_addfast_cc
        {
            iof(output,4*i+j)=cos_vals[j]*iof(input,4*i+j)-sin_vals[j]*qof(input,4*i+j);
            qof(output,4*i+j)=sin_vals[j]*iof(input,4*i+j)+cos_vals[j]*qof(input,4*i+j);
        }
        cos_start = cos_vals[3];
        sin_start = sin_vals[3];
    }
    starting_phase+=input_size*d->phase_increment;
    while(starting_phase>PI) starting_phase-=2*PI;
    while(starting_phase<-PI) starting_phase+=2*PI;
    return starting_phase;
}
#endif

#endif

#ifdef NEON_OPTS
#pragma message "Manual NEON optimizations are ON: we have a faster fir_decimate_cc now."

//max help: http://community.arm.com/groups/android-community/blog/2015/03/27/arm-neon-programming-quick-reference

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
        register float* pinput=(float*)&(input[i]);
        register float* ptaps=taps;
        register float* ptaps_end=taps+taps_length;
        float quad_acciq [8];


/*
q0, q1: input signal I sample and Q sample
q2:     taps
q4, q5: accumulator for I branch and Q branch (will be the output)
*/

        asm volatile(
            "       veor q4, q4\n\t"
            "       veor q5, q5\n\t"
            "for_fdccasm: vld2.32   {q0-q1}, [%[pinput]]!\n\t" //load q0 and q1 directly from the memory address stored in pinput, with interleaving (so that we get the I samples in q0 and the Q samples in q1), also increment the memory address in pinput (hence the "!" mark) //http://community.arm.com/groups/processors/blog/2010/03/17/coding-for-neon--part-1-load-and-stores
            "       vld1.32 {q2}, [%[ptaps]]!\n\t"
            "       vmla.f32 q4, q0, q2\n\t" //quad_acc_i += quad_input_i * quad_taps_1 //http://stackoverflow.com/questions/3240440/how-to-use-the-multiply-and-accumulate-intrinsics-in-arm-cortex-a8 //http://infocenter.arm.com/help/index.jsp?topic=/com.arm.doc.dui0489e/CIHEJBIE.html
            "       vmla.f32 q5, q1, q2\n\t" //quad_acc_q += quad_input_q * quad_taps_1
            "       cmp %[ptaps], %[ptaps_end]\n\t" //if(ptaps != ptaps_end)
            "       bcc for_fdccasm\n\t"            //  then goto for_fdcasm
            "       vst1.32 {q4}, [%[quad_acci]]\n\t" //if the loop is finished, store the two accumulators in memory
            "       vst1.32 {q5}, [%[quad_accq]]\n\t"
        :
            [pinput]"+r"(pinput), [ptaps]"+r"(ptaps) //output operand list
        :
            [ptaps_end]"r"(ptaps_end), [quad_acci]"r"(quad_acciq), [quad_accq]"r"(quad_acciq+4) //input operand list
        :
            "memory", "q0", "q1", "q2", "q4", "q5", "cc" //clobber list
        );
        //original for loops for reference:
        //for(int ti=0; ti<taps_length; ti++) acci += (iof(input,i+ti)) * taps[ti]; //@fir_decimate_cc: i loop
        //for(int ti=0; ti<taps_length; ti++) accq += (qof(input,i+ti)) * taps[ti]; //@fir_decimate_cc: q loop

        //for(int n=0;n<8;n++) fprintf(stderr, "\n>> [%d] %g \n", n, quad_acciq[n]);
        iof(output,oi)=quad_acciq[0]+quad_acciq[1]+quad_acciq[2]+quad_acciq[3]; //we're still not ready, as we have to add up the contents of a quad accumulator register to get a single accumulated value
        qof(output,oi)=quad_acciq[4]+quad_acciq[5]+quad_acciq[6]+quad_acciq[7];
        oi++;
    }
    return oi;
}

#else

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

#endif

/*
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
        int taps_halflength = taps_length/2;
        for(int ti=0; ti<taps_halflength; ti++) acci += (iof(input,i+ti)+iof(input,i+taps_length-ti)) * taps[ti]; //@fir_decimate_cc: i loop
        float accq=0;
        for(int ti=0; ti<taps_halflength; ti++) accq += (qof(input,i+ti)+qof(input,i+taps_length-ti)) * taps[ti]; //@fir_decimate_cc: q loop
        iof(output,oi)=acci+iof(input,i+taps_halflength)*taps[taps_halflength];
        qof(output,oi)=accq+qof(input,i+taps_halflength)*taps[taps_halflength];
        oi++;
    }
    return oi;
}
*/

int fir_interpolate_cc(complexf *input, complexf *output, int input_size, int interpolation, float *taps, int taps_length)
{
    //i:  input index
    //oi: output index
    //ti: tap index
    //ti: secondary index (inside filter function)
    //ip: interpolation phase (0 <= ip < interpolation)
    int oi=0;
    for(int i=0; i<input_size; i++) //@fir_interpolate_cc: outer loop
    {
        if(i*interpolation + (interpolation-1) + taps_length > input_size*interpolation) break;
        for(int ip=0; ip<interpolation; ip++)
        {
            float acci=0;
            float accq=0;
            //int tistart = (interpolation-ip)%interpolation; 
            int tistart = (interpolation-ip); //why does this work? why don't we need the % part?
            for(int ti=tistart, si=0; ti<taps_length; (ti+=interpolation), (si++)) acci += (iof(input,i+si)) * taps[ti]; //@fir_interpolate_cc: i loop
            for(int ti=tistart, si=0; ti<taps_length; (ti+=interpolation), (si++)) accq += (qof(input,i+si)) * taps[ti]; //@fir_interpolate_cc: q loop
            iof(output,oi)=acci;
            qof(output,oi)=accq;
            oi++;
        }
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
        for(int i=0; i<(taps_length-delayi)/interpolation; i++) //@rational_resampler_ff (inner loop)
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

old_fractional_decimator_ff_t old_fractional_decimator_ff(float* input, float* output, int input_size, float rate, float *taps, int taps_length, old_fractional_decimator_ff_t d)
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

fractional_decimator_ff_t fractional_decimator_ff_init(float rate, int num_poly_points, float* taps, int taps_length)
{
    fractional_decimator_ff_t d;
    d.num_poly_points = num_poly_points&~1; //num_poly_points needs to be even!
    d.poly_precalc_denomiator = (float*)malloc(d.num_poly_points*sizeof(float));
    //x0..x3
    //-1,0,1,2
    //-(4/2)+1
    //x0..x5
    //-2,-1,0,1,2,3
    d.xifirst=-(num_poly_points/2)+1, d.xilast=num_poly_points/2;
    int id = 0; //index in poly_precalc_denomiator
    for(int xi=d.xifirst;xi<=d.xilast;xi++)
    {
        d.poly_precalc_denomiator[id]=1;
        for(int xj=d.xifirst;xj<=d.xilast;xj++)
        {
            if(xi!=xj) d.poly_precalc_denomiator[id] *= (xi-xj); //poly_precalc_denomiator could be integer as well. But that would later add a necessary conversion.
        }
        id++;
    }
    d.where=-d.xifirst;
    d.coeffs_buf=(float*)malloc(d.num_poly_points*sizeof(float)); 
    d.filtered_buf=(float*)malloc(d.num_poly_points*sizeof(float)); 
    //d.last_inputs_circbuf = (float)malloc(d.num_poly_points*sizeof(float));
    //d.last_inputs_startsat = 0; 
    //d.last_inputs_samplewhere = -1;
    //for(int i=0;i<num_poly_points; i++) d.last_inputs_circbuf[i] = 0;
    d.rate = rate;
    d.taps = taps;
    d.taps_length = taps_length;
    d.input_processed = 0;
    return d;
}

#define DEBUG_ASSERT 1
void fractional_decimator_ff(float* input, float* output, int input_size, fractional_decimator_ff_t* d)
{
    //This routine can handle floating point decimation rates.
    //It applies polynomial interpolation to samples that are taken into consideration from a pre-filtered input.
    //The pre-filter can be switched off by applying taps=NULL.
    //fprintf(stderr, "drate=%f\n", d->rate);
    if(DEBUG_ASSERT) assert(d->rate > 1.0); 
    if(DEBUG_ASSERT) assert(d->where >= -d->xifirst);
    int oi=0; //output index
    int index_high; 
#define FD_INDEX_LOW (index_high-1)
    //we optimize to calculate ceilf(where) only once every iteration, so we do it here:
    for(;(index_high=ceilf(d->where))+d->num_poly_points+d->taps_length<input_size;d->where+=d->rate) //@fractional_decimator_ff
    {
        //d->num_poly_points above is theoretically more than we could have here, but this makes the spectrum look good
        int sxifirst = FD_INDEX_LOW + d->xifirst; 
        int sxilast = FD_INDEX_LOW + d->xilast; 
        if(d->taps) 
            for(int wi=0;wi<d->num_poly_points;wi++) d->filtered_buf[wi] = fir_one_pass_ff(input+FD_INDEX_LOW+wi, d->taps, d->taps_length);
        else
            for(int wi=0;wi<d->num_poly_points;wi++) d->filtered_buf[wi] = *(input+FD_INDEX_LOW+wi);
        int id=0;
        float xwhere = d->where - FD_INDEX_LOW;
        for(int xi=d->xifirst;xi<=d->xilast;xi++)
        {
            d->coeffs_buf[id]=1;
            for(int xj=d->xifirst;xj<=d->xilast;xj++)
            {
                if(xi!=xj) d->coeffs_buf[id] *= (xwhere-xj);
            }
            id++;       
        }
        float acc = 0;
        for(int i=0;i<d->num_poly_points;i++)
        {
            acc += (d->coeffs_buf[i]/d->poly_precalc_denomiator[i])*d->filtered_buf[i];  //(xnom/xden)*yn
        }
        output[oi++]=acc;
    }
    d->input_processed = FD_INDEX_LOW + d->xifirst;
    d->where -= d->input_processed;
    d->output_size = oi;
}

/*
 * Some notes to myself on the circular buffer I wanted to implement here:
        int last_input_samplewhere_shouldbe = (index_high-1)+xifirst;
        int last_input_offset = last_input_samplewhere_shouldbe - d->last_input_samplewhere;
        if(last_input_offset < num_poly_points)
        {
            //if we can move the last_input circular buffer, we move, and add the new samples at the end
            d->last_inputs_startsat += last_input_offset;
            d->last_inputs_startsat %= num_poly_points;
            int num_copied_samples = 0;
            for(int i=0; i<last_input_offset; i++)
            {
                d->last_inputs_circbuf[i]=
            }
            d->last_input_samplewhere = d->las
        }
    However, I think I should just rather do a continuous big buffer.
*/

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

//#define FASTAGC_MAX_GAIN (65e3)
#define FASTAGC_MAX_GAIN 50

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
    //fprintf(stderr, "target_gain: %g\n",target_gain);

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
        output[i]=(temp[i])?fmdemod_quadri_K*output[i]/temp[i]:0;
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

float get_power_f(float* input, int input_size, int decimation)
{
    float acc = 0;
    for(int i=0;i<input_size;i+=decimation)
    {
        acc += (input[i]*input[i])/input_size;
    }
    return acc;
}

float get_power_c(complexf* input, int input_size, int decimation)
{
    float acc = 0;
    for(int i=0;i<input_size;i+=decimation)
    {
        acc += (iof(input,i)*iof(input,i)+qof(input,i)*qof(input,i))/input_size;
    }
    return acc;
}

/*
  __  __           _       _       _ 
 |  \/  |         | |     | |     |  |
 | \  / | ___   __| |_   _| | __ _|  |_ ___  _ __ ___
 | |\/| |/ _ \ / _` | | | | |/ _` |  __/ _ \| '__/ __|
 | |  | | (_) | (_| | |_| | | (_| |  || (_) | |  \__ \
 |_|  |_|\___/ \__,_|\__,_|_|\__,_|\ __\___/|_|  |___/

*/

void add_dcoffset_cc(complexf* input, complexf* output, int input_size)
{
    for(int i=0;i<input_size;i++) iof(output,i)=0.5+iof(input,i)/2;
    for(int i=0;i<input_size;i++) qof(output,i)=qof(input,i)/2;
}

float fmmod_fc(float* input, complexf* output, int input_size, float last_phase)
{
    float phase=last_phase;
    for(int i=0;i<input_size;i++)
    {
        phase+=input[i]*PI;
        while(phase>PI) phase-=2*PI;
        while(phase<=-PI) phase+=2*PI;
        iof(output,i)=cos(phase);
        qof(output,i)=sin(phase);
    }
    return phase;
}

void fixed_amplitude_cc(complexf* input, complexf* output, int input_size, float new_amplitude)
{
    for(int i=0;i<input_size;i++)
    {
        //float phase=atan2(iof(input,i),qof(input,i));
        //iof(output,i)=cos(phase)*amp;
        //qof(output,i)=sin(phase)*amp;

        //A faster solution:
        float amplitude_now = sqrt(iof(input,i)*iof(input,i)+qof(input,i)*qof(input,i));
        float gain = (amplitude_now > 0) ? new_amplitude / amplitude_now : 0;
        iof(output,i)=iof(input,i)*gain;
        qof(output,i)=qof(input,i)*gain;
    }
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

float *precalculate_window(int size, window_t window)
{
    float (*window_function)(float)=firdes_get_window_kernel(window);
    float *windowt;
    windowt = malloc(sizeof(float) * size);
    for(int i=0;i<size;i++) //@precalculate_window
    {
        float rate=(float)i/(size-1);
        windowt[i] = window_function(2.0*rate+1.0);
    }
    return windowt;
}

void apply_precalculated_window_c(complexf* input, complexf* output, int size, float *windowt)
{
    for(int i=0;i<size;i++) //@apply_precalculated_window_c
    {
        iof(output,i)=iof(input,i)*windowt[i];
        qof(output,i)=qof(input,i)*windowt[i];
    }
}

void apply_precalculated_window_f(float* input, float* output, int size, float *windowt)
{
	for(int i=0;i<size;i++) //@apply_precalculated_window_f
	{
		output[i] = input[i] * windowt[i];
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

void accumulate_power_cf(complexf* input, float* output, int size)
{
    for(int i=0;i<size;i++) output[i] += iof(input,i)*iof(input,i) + qof(input,i)*qof(input,i); //@logpower_cf: pass 1
}

void log_ff(float* input, float* output, int size, float add_db) {
    for(int i=0;i<size;i++) output[i]=log10(input[i]); //@logpower_cf: pass 2

    for(int i=0;i<size;i++) output[i]=10*output[i]+add_db; //@logpower_cf: pass 3
}

float total_logpower_cf(complexf* input, int input_size)
{
    float acc = 0; 
    for(int i=0;i<input_size;i++) acc+=(iof(input,i)*iof(input,i) + qof(input,i)*qof(input,i));
    return 10*log10(acc/input_size);
}

/*
  _____  _       _ _        _       _                          _
 |  __ \(_)     (_) |      | |     | |                        | |
 | |  | |_  __ _ _| |_ __ _| |   __| | ___ _ __ ___   ___   __| |
 | |  | | |/ _` | | __/ _` | |  / _` |/ _ \ '_ ` _ \ / _ \ / _` |
 | |__| | | (_| | | || (_| | | | (_| |  __/ | | | | | (_) | (_| |
 |_____/|_|\__, |_|\__\__,_|_|  \__,_|\___|_| |_| |_|\___/ \__,_|
            __/ |
           |___/
*/

psk31_varicode_item_t psk31_varicode_items[] =
{
    { .code = 0b1010101011, .bitcount=10,   .ascii=0x00 }, //NUL, null
    { .code = 0b1011011011, .bitcount=10,   .ascii=0x01 }, //SOH, start of heading
    { .code = 0b1011101101, .bitcount=10,   .ascii=0x02 }, //STX, start of text
    { .code = 0b1101110111, .bitcount=10,   .ascii=0x03 }, //ETX, end of text
    { .code = 0b1011101011, .bitcount=10,   .ascii=0x04 }, //EOT, end of transmission
    { .code = 0b1101011111, .bitcount=10,   .ascii=0x05 }, //ENQ, enquiry
    { .code = 0b1011101111, .bitcount=10,   .ascii=0x06 }, //ACK, acknowledge
    { .code = 0b1011111101, .bitcount=10,   .ascii=0x07 }, //BEL, bell
    { .code = 0b1011111111, .bitcount=10,   .ascii=0x08 }, //BS, backspace
    { .code = 0b11101111,   .bitcount=8,    .ascii=0x09 }, //TAB, horizontal tab
    { .code = 0b11101,      .bitcount=5,    .ascii=0x0a }, //LF, NL line feed, new line
    { .code = 0b1101101111, .bitcount=10,   .ascii=0x0b }, //VT, vertical tab
    { .code = 0b1011011101, .bitcount=10,   .ascii=0x0c }, //FF, NP form feed, new page
    { .code = 0b11111,      .bitcount=5,    .ascii=0x0d }, //CR, carriage return (overwrite)
    { .code = 0b1101110101, .bitcount=10,   .ascii=0x0e }, //SO, shift out
    { .code = 0b1110101011, .bitcount=10,   .ascii=0x0f }, //SI, shift in
    { .code = 0b1011110111, .bitcount=10,   .ascii=0x10 }, //DLE, data link escape
    { .code = 0b1011110101, .bitcount=10,   .ascii=0x11 }, //DC1, device control 1
    { .code = 0b1110101101, .bitcount=10,   .ascii=0x12 }, //DC2, device control 2
    { .code = 0b1110101111, .bitcount=10,   .ascii=0x13 }, //DC3, device control 3
    { .code = 0b1101011011, .bitcount=10,   .ascii=0x14 }, //DC4, device control 4
    { .code = 0b1101101011, .bitcount=10,   .ascii=0x15 }, //NAK, negative acknowledge
    { .code = 0b1101101101, .bitcount=10,   .ascii=0x16 }, //SYN, synchronous idle
    { .code = 0b1101010111, .bitcount=10,   .ascii=0x17 }, //ETB, end of trans. block
    { .code = 0b1101111011, .bitcount=10,   .ascii=0x18 }, //CAN, cancel
    { .code = 0b1101111101, .bitcount=10,   .ascii=0x19 }, //EM, end of medium
    { .code = 0b1110110111, .bitcount=10,   .ascii=0x1a }, //SUB, substitute
    { .code = 0b1101010101, .bitcount=10,   .ascii=0x1b }, //ESC, escape
    { .code = 0b1101011101, .bitcount=10,   .ascii=0x1c }, //FS, file separator
    { .code = 0b1110111011, .bitcount=10,   .ascii=0x1d }, //GS, group separator
    { .code = 0b1011111011, .bitcount=10,   .ascii=0x1e }, //RS, record separator
    { .code = 0b1101111111, .bitcount=10,   .ascii=0x1f }, //US, unit separator
    { .code = 0b1,          .bitcount=1,    .ascii=0x20 }, //szóköz
    { .code = 0b111111111,  .bitcount=9,    .ascii=0x21 }, //!
    { .code = 0b101011111,  .bitcount=9,    .ascii=0x22 }, //"
    { .code = 0b111110101,  .bitcount=9,    .ascii=0x23 }, //#
    { .code = 0b111011011,  .bitcount=9,    .ascii=0x24 }, //$
    { .code = 0b1011010101, .bitcount=10,   .ascii=0x25 }, //%
    { .code = 0b1010111011, .bitcount=10,   .ascii=0x26 }, //&
    { .code = 0b101111111,  .bitcount=9,    .ascii=0x27 }, //'
    { .code = 0b11111011,   .bitcount=8,    .ascii=0x28 }, //(
    { .code = 0b11110111,   .bitcount=8,    .ascii=0x29 }, //)
    { .code = 0b101101111,  .bitcount=9,    .ascii=0x2a }, //*
    { .code = 0b111011111,  .bitcount=9,    .ascii=0x2b }, //+
    { .code = 0b1110101,    .bitcount=7,    .ascii=0x2c }, //,
    { .code = 0b110101,     .bitcount=6,    .ascii=0x2d }, //-
    { .code = 0b1010111,    .bitcount=7,    .ascii=0x2e }, //.
    { .code = 0b110101111,  .bitcount=9,    .ascii=0x2f }, ///
    { .code = 0b10110111,   .bitcount=8,    .ascii=0x30 }, //0
    { .code = 0b10111101,   .bitcount=8,    .ascii=0x31 }, //1
    { .code = 0b11101101,   .bitcount=8,    .ascii=0x32 }, //2
    { .code = 0b11111111,   .bitcount=8,    .ascii=0x33 }, //3
    { .code = 0b101110111,  .bitcount=9,    .ascii=0x34 }, //4
    { .code = 0b101011011,  .bitcount=9,    .ascii=0x35 }, //5
    { .code = 0b101101011,  .bitcount=9,    .ascii=0x36 }, //6
    { .code = 0b110101101,  .bitcount=9,    .ascii=0x37 }, //7
    { .code = 0b110101011,  .bitcount=9,    .ascii=0x38 }, //8
    { .code = 0b110110111,  .bitcount=9,    .ascii=0x39 }, //9
    { .code = 0b11110101,   .bitcount=8,    .ascii=0x3a }, //:
    { .code = 0b110111101,  .bitcount=9,    .ascii=0x3b }, //;
    { .code = 0b111101101,  .bitcount=9,    .ascii=0x3c }, //<
    { .code = 0b1010101,    .bitcount=7,    .ascii=0x3d }, //=
    { .code = 0b111010111,  .bitcount=9,    .ascii=0x3e }, //>
    { .code = 0b1010101111, .bitcount=10,   .ascii=0x3f }, //?
    { .code = 0b1010111101, .bitcount=10,   .ascii=0x40 }, //@
    { .code = 0b1111101,    .bitcount=7,    .ascii=0x41 }, //A
    { .code = 0b11101011,   .bitcount=8,    .ascii=0x42 }, //B
    { .code = 0b10101101,   .bitcount=8,    .ascii=0x43 }, //C
    { .code = 0b10110101,   .bitcount=8,    .ascii=0x44 }, //D
    { .code = 0b1110111,    .bitcount=7,    .ascii=0x45 }, //E
    { .code = 0b11011011,   .bitcount=8,    .ascii=0x46 }, //F
    { .code = 0b11111101,   .bitcount=8,    .ascii=0x47 }, //G
    { .code = 0b101010101,  .bitcount=9,    .ascii=0x48 }, //H
    { .code = 0b1111111,    .bitcount=7,    .ascii=0x49 }, //I
    { .code = 0b111111101,  .bitcount=9,    .ascii=0x4a }, //J
    { .code = 0b101111101,  .bitcount=9,    .ascii=0x4b }, //K
    { .code = 0b11010111,   .bitcount=8,    .ascii=0x4c }, //L
    { .code = 0b10111011,   .bitcount=8,    .ascii=0x4d }, //M
    { .code = 0b11011101,   .bitcount=8,    .ascii=0x4e }, //N
    { .code = 0b10101011,   .bitcount=8,    .ascii=0x4f }, //O
    { .code = 0b11010101,   .bitcount=8,    .ascii=0x50 }, //P
    { .code = 0b111011101,  .bitcount=9,    .ascii=0x51 }, //Q
    { .code = 0b10101111,   .bitcount=8,    .ascii=0x52 }, //R
    { .code = 0b1101111,    .bitcount=7,    .ascii=0x53 }, //S
    { .code = 0b1101101,    .bitcount=7,    .ascii=0x54 }, //T
    { .code = 0b101010111,  .bitcount=9,    .ascii=0x55 }, //U
    { .code = 0b110110101,  .bitcount=9,    .ascii=0x56 }, //V
    { .code = 0b101011101,  .bitcount=9,    .ascii=0x57 }, //W
    { .code = 0b101110101,  .bitcount=9,    .ascii=0x58 }, //X
    { .code = 0b101111011,  .bitcount=9,    .ascii=0x59 }, //Y
    { .code = 0b1010101101, .bitcount=10,   .ascii=0x5a }, //Z
    { .code = 0b111110111,  .bitcount=9,    .ascii=0x5b }, //[
    { .code = 0b111101111,  .bitcount=9,    .ascii=0x5c }, //backslash
    { .code = 0b111111011,  .bitcount=9,    .ascii=0x5d }, //]
    { .code = 0b1010111111, .bitcount=10,   .ascii=0x5e }, //^
    { .code = 0b101101101,  .bitcount=9,    .ascii=0x5f }, //_
    { .code = 0b1011011111, .bitcount=10,   .ascii=0x60 }, //`
    { .code = 0b1011,       .bitcount=4,    .ascii=0x61 }, //a
    { .code = 0b1011111,    .bitcount=7,    .ascii=0x62 }, //b
    { .code = 0b101111,     .bitcount=6,    .ascii=0x63 }, //c
    { .code = 0b101101,     .bitcount=6,    .ascii=0x64 }, //d
    { .code = 0b11,         .bitcount=2,    .ascii=0x65 }, //e
    { .code = 0b111101,     .bitcount=6,    .ascii=0x66 }, //f
    { .code = 0b1011011,    .bitcount=7,    .ascii=0x67 }, //g
    { .code = 0b101011,     .bitcount=6,    .ascii=0x68 }, //h
    { .code = 0b1101,       .bitcount=4,    .ascii=0x69 }, //i
    { .code = 0b111101011,  .bitcount=9,    .ascii=0x6a }, //j
    { .code = 0b10111111,   .bitcount=8,    .ascii=0x6b }, //k
    { .code = 0b11011,      .bitcount=5,    .ascii=0x6c }, //l
    { .code = 0b111011,     .bitcount=6,    .ascii=0x6d }, //m
    { .code = 0b1111,       .bitcount=4,    .ascii=0x6e }, //n
    { .code = 0b111,        .bitcount=3,    .ascii=0x6f }, //o
    { .code = 0b111111,     .bitcount=6,    .ascii=0x70 }, //p
    { .code = 0b110111111,  .bitcount=9,    .ascii=0x71 }, //q
    { .code = 0b10101,      .bitcount=5,    .ascii=0x72 }, //r
    { .code = 0b10111,      .bitcount=5,    .ascii=0x73 }, //s
    { .code = 0b101,        .bitcount=3,    .ascii=0x74 }, //t
    { .code = 0b110111,     .bitcount=6,    .ascii=0x75 }, //u
    { .code = 0b1111011,    .bitcount=7,    .ascii=0x76 }, //v
    { .code = 0b1101011,    .bitcount=7,    .ascii=0x77 }, //w
    { .code = 0b11011111,   .bitcount=8,    .ascii=0x78 }, //x
    { .code = 0b1011101,    .bitcount=7,    .ascii=0x79 }, //y
    { .code = 0b111010101,  .bitcount=9,    .ascii=0x7a }, //z
    { .code = 0b1010110111, .bitcount=10,   .ascii=0x7b }, //{
    { .code = 0b110111011,  .bitcount=9,    .ascii=0x7c }, //|
    { .code = 0b1010110101, .bitcount=10,   .ascii=0x7d }, //}
    { .code = 0b1011010111, .bitcount=10,   .ascii=0x7e }, //~
    { .code = 0b1110110101, .bitcount=10,   .ascii=0x7f }, //DEL
};

unsigned long long psk31_varicode_masklen_helper[] =
{
    0b0000000000000000000000000000000000000000000000000000000000000000,
    0b0000000000000000000000000000000000000000000000000000000000000001,
    0b0000000000000000000000000000000000000000000000000000000000000011,
    0b0000000000000000000000000000000000000000000000000000000000000111,
    0b0000000000000000000000000000000000000000000000000000000000001111,
    0b0000000000000000000000000000000000000000000000000000000000011111,
    0b0000000000000000000000000000000000000000000000000000000000111111,
    0b0000000000000000000000000000000000000000000000000000000001111111,
    0b0000000000000000000000000000000000000000000000000000000011111111,
    0b0000000000000000000000000000000000000000000000000000000111111111,
    0b0000000000000000000000000000000000000000000000000000001111111111,
    0b0000000000000000000000000000000000000000000000000000011111111111,
    0b0000000000000000000000000000000000000000000000000000111111111111,
    0b0000000000000000000000000000000000000000000000000001111111111111,
    0b0000000000000000000000000000000000000000000000000011111111111111,
    0b0000000000000000000000000000000000000000000000000111111111111111,
    0b0000000000000000000000000000000000000000000000001111111111111111,
    0b0000000000000000000000000000000000000000000000011111111111111111,
    0b0000000000000000000000000000000000000000000000111111111111111111,
    0b0000000000000000000000000000000000000000000001111111111111111111,
    0b0000000000000000000000000000000000000000000011111111111111111111,
    0b0000000000000000000000000000000000000000000111111111111111111111,
    0b0000000000000000000000000000000000000000001111111111111111111111,
    0b0000000000000000000000000000000000000000011111111111111111111111,
    0b0000000000000000000000000000000000000000111111111111111111111111,
    0b0000000000000000000000000000000000000001111111111111111111111111,
    0b0000000000000000000000000000000000000011111111111111111111111111,
    0b0000000000000000000000000000000000000111111111111111111111111111,
    0b0000000000000000000000000000000000001111111111111111111111111111,
    0b0000000000000000000000000000000000011111111111111111111111111111,
    0b0000000000000000000000000000000000111111111111111111111111111111,
    0b0000000000000000000000000000000001111111111111111111111111111111,
    0b0000000000000000000000000000000011111111111111111111111111111111,
    0b0000000000000000000000000000000111111111111111111111111111111111,
    0b0000000000000000000000000000001111111111111111111111111111111111,
    0b0000000000000000000000000000011111111111111111111111111111111111,
    0b0000000000000000000000000000111111111111111111111111111111111111,
    0b0000000000000000000000000001111111111111111111111111111111111111,
    0b0000000000000000000000000011111111111111111111111111111111111111,
    0b0000000000000000000000000111111111111111111111111111111111111111,
    0b0000000000000000000000001111111111111111111111111111111111111111,
    0b0000000000000000000000011111111111111111111111111111111111111111,
    0b0000000000000000000000111111111111111111111111111111111111111111,
    0b0000000000000000000001111111111111111111111111111111111111111111,
    0b0000000000000000000011111111111111111111111111111111111111111111,
    0b0000000000000000000111111111111111111111111111111111111111111111,
    0b0000000000000000001111111111111111111111111111111111111111111111,
    0b0000000000000000011111111111111111111111111111111111111111111111,
    0b0000000000000000111111111111111111111111111111111111111111111111,
    0b0000000000000001111111111111111111111111111111111111111111111111,
    0b0000000000000011111111111111111111111111111111111111111111111111,
    0b0000000000000111111111111111111111111111111111111111111111111111,
    0b0000000000001111111111111111111111111111111111111111111111111111,
    0b0000000000011111111111111111111111111111111111111111111111111111,
    0b0000000000111111111111111111111111111111111111111111111111111111,
    0b0000000001111111111111111111111111111111111111111111111111111111,
    0b0000000011111111111111111111111111111111111111111111111111111111,
    0b0000000111111111111111111111111111111111111111111111111111111111,
    0b0000001111111111111111111111111111111111111111111111111111111111,
    0b0000011111111111111111111111111111111111111111111111111111111111,
    0b0000111111111111111111111111111111111111111111111111111111111111,
    0b0001111111111111111111111111111111111111111111111111111111111111,
    0b0011111111111111111111111111111111111111111111111111111111111111,
    0b0111111111111111111111111111111111111111111111111111111111111111
};

const int n_psk31_varicode_items = sizeof(psk31_varicode_items) / sizeof(psk31_varicode_item_t);

char psk31_varicode_decoder_push(unsigned long long* status_shr, unsigned char symbol)
{
    *status_shr=((*status_shr)<<1)|(!!symbol); //shift new bit in shift register
    //fprintf(stderr,"*status_shr = %llx\n", *status_shr);
    if((*status_shr)&0xFFF==0) return 0;
    for(int i=0;i<n_psk31_varicode_items;i++)
    {
        //fprintf(stderr,"| i = %d | %llx ?= %llx | bitsall = %d\n", i, psk31_varicode_items[i].code<<2, (*status_shr)&psk31_varicode_masklen_helper[(psk31_varicode_items[i].bitcount+4)&63], (psk31_varicode_items[i].bitcount+4)&63);
        if((psk31_varicode_items[i].code<<2)==((*status_shr)&psk31_varicode_masklen_helper[(psk31_varicode_items[i].bitcount+4)&63]))
            { /*fprintf(stderr,">>>>>>>>> %d %x %c\n", i,  psk31_varicode_items[i].ascii,  psk31_varicode_items[i].ascii);*/ return psk31_varicode_items[i].ascii; }

    }
    return 0;
}

void psk31_varicode_encoder_u8_u8(unsigned char* input, unsigned char* output, int input_size, int output_max_size, int* input_processed, int* output_size)
{
    (*output_size)=0;
    for((*input_processed)=0; (*input_processed)<input_size; (*input_processed)++)
    {
        //fprintf(stderr, "ii = %d, input_size = %d, output_max_size = %d\n", *input_processed, input_size, output_max_size);
        for(int ci=0; ci<n_psk31_varicode_items; ci++) //ci: character index
        {
            psk31_varicode_item_t current_varicode = psk31_varicode_items[ci];
            if(input[*input_processed]==current_varicode.ascii)
            {
                //fprintf(stderr, "ci = %d\n", ci);
                if(output_max_size<current_varicode.bitcount+2) return;
                for(int bi=0; bi<current_varicode.bitcount+2; bi++) //bi: bit index
                {
                    //fprintf(stderr, "bi = %d\n", bi);
                    output[*output_size] = (bi<current_varicode.bitcount) ? (psk31_varicode_items[ci].code>>(current_varicode.bitcount-bi-1))&1 : 0;
                    (*output_size)++;
                    output_max_size--;
                }
                break;
            }
        }
    }
}

rtty_baudot_item_t rtty_baudot_items[] =
{
    { .code = 0b00000, .ascii_letter=0,     .ascii_figure=0 },
    { .code = 0b10000, .ascii_letter='E',   .ascii_figure='3' },
    { .code = 0b01000, .ascii_letter='\n',  .ascii_figure='\n' },
    { .code = 0b11000, .ascii_letter='A',   .ascii_figure='-' },
    { .code = 0b00100, .ascii_letter=' ',   .ascii_figure=' ' },
    { .code = 0b10100, .ascii_letter='S',   .ascii_figure='\'' },
    { .code = 0b01100, .ascii_letter='I',   .ascii_figure='8' },
    { .code = 0b11100, .ascii_letter='U',   .ascii_figure='7' },
    { .code = 0b00010, .ascii_letter='\r',  .ascii_figure='\r' },
    { .code = 0b10010, .ascii_letter='D',   .ascii_figure='#' },
    { .code = 0b01010, .ascii_letter='R',   .ascii_figure='4' },
    { .code = 0b11010, .ascii_letter='J',   .ascii_figure='\a' },
    { .code = 0b00110, .ascii_letter='N',   .ascii_figure=',' },
    { .code = 0b10110, .ascii_letter='F',   .ascii_figure='@' },
    { .code = 0b01110, .ascii_letter='C',   .ascii_figure=':' },
    { .code = 0b11110, .ascii_letter='K',   .ascii_figure='(' },
    { .code = 0b00001, .ascii_letter='T',   .ascii_figure='5' },
    { .code = 0b10001, .ascii_letter='Z',   .ascii_figure='+' },
    { .code = 0b01001, .ascii_letter='L',   .ascii_figure=')' },
    { .code = 0b11001, .ascii_letter='W',   .ascii_figure='2' },
    { .code = 0b00101, .ascii_letter='H',   .ascii_figure='$' },
    { .code = 0b10101, .ascii_letter='Y',   .ascii_figure='6' },
    { .code = 0b01101, .ascii_letter='P',   .ascii_figure='0' },
    { .code = 0b11101, .ascii_letter='Q',   .ascii_figure='1' },
    { .code = 0b00011, .ascii_letter='O',   .ascii_figure='9' },
    { .code = 0b10011, .ascii_letter='B',   .ascii_figure='?' },
    { .code = 0b01011, .ascii_letter='G',   .ascii_figure='*' },
    { .code = 0b00111, .ascii_letter='M',   .ascii_figure='.' },
    { .code = 0b10111, .ascii_letter='X',   .ascii_figure='/' },
    { .code = 0b01111, .ascii_letter='V',   .ascii_figure='=' }
};

const int n_rtty_baudot_items = sizeof(rtty_baudot_items) / sizeof(rtty_baudot_item_t);

char rtty_baudot_decoder_lookup(unsigned char* fig_mode, unsigned char c)
{
    if(c==RTTY_FIGURE_MODE_SELECT_CODE) { *fig_mode=1; return 0; }
    if(c==RTTY_LETTER_MODE_SELECT_CODE) { *fig_mode=0; return 0; }
    for(int i=0;i<n_rtty_baudot_items;i++)
        if(rtty_baudot_items[i].code==c)
            return (*fig_mode) ? rtty_baudot_items[i].ascii_figure : rtty_baudot_items[i].ascii_letter;
    return 0;
}

char rtty_baudot_decoder_push(rtty_baudot_decoder_t* s, unsigned char symbol)
{
    //For RTTY waveforms, check this: http://www.ham.hu/radiosatvitel/szoveg/RTTY/kepek/rtty.gif
    //RTTY is much like an UART data transfer with 1 start bit, 5 data bits and 1 stop bit.
    //The start pulse and stop pulse are used for synchronization.
    symbol=!!symbol; //We want symbol to be 0 or 1.
    switch(s->state)
    {
    case RTTY_BAUDOT_WAITING_STOP_PULSE:
        if(symbol==1) { s->state = RTTY_BAUDOT_WAITING_START_PULSE; if(s->character_received) return rtty_baudot_decoder_lookup(&s->fig_mode, s->shr&31); }
            //If the character data is followed by a stop pulse, then we go on to wait for the next character.
        else  s->character_received = 0;
            //The character should be followed by a stop pulse. If the stop pulse is missing, that is certainly an error.
            //In that case, we remove forget the character we just received.
        break;
    case RTTY_BAUDOT_WAITING_START_PULSE:
        s->character_received = 0;
        if(symbol==0) { s->state = RTTY_BAUDOT_RECEIVING_DATA; s->shr = s->bit_cntr = 0; }
            //Any number of high bits can come after each other, until interrupted with a low bit (start pulse) to indicate
            //the beginning of a new character. If we get this start pulse, we go on to wait for the characters. We also
            //clear the variables used for counting (bit_cntr) and storing (shr) the data bits.
        break;
    case RTTY_BAUDOT_RECEIVING_DATA:
        s->shr = (s->shr<<1)|(!!symbol);
            //We store 5 bits into our shift register
        if(s->bit_cntr++==4) { s->state = RTTY_BAUDOT_WAITING_STOP_PULSE; s->character_received = 1; }
            //If this is the 5th bit stored, then we wait for the stop pulse.
        break;
    default: break;
    }
    return 0;
}

#define DEBUG_SERIAL_LINE_DECODER 0

//What has not been checked:
//  behaviour on 1.5 stop bits
//  check all exit conditions

void serial_line_decoder_f_u8(serial_line_t* s, float* input, unsigned char* output, int input_size)
{
    static int abs_samples_helper = 0;
    abs_samples_helper += s->input_used;
    int iabs_samples_helper = abs_samples_helper;
    s->output_size = 0;
    s->input_used = 0;
    short* output_s = (short*)output;
    unsigned* output_u = (unsigned*)output;
    for(;;)
    {
        //we find the start bit (first negative edge on the line)
        int startbit_start = -1;
        int i;
        for(i=1;i<input_size;i++) if(input[i] < 0 && input[i-1] > 0) { startbit_start=i; break; }

        if(startbit_start == -1) { s->input_used += i; DEBUG_SERIAL_LINE_DECODER && fprintf(stderr,"sld:startbit_not_found (+%d)\n", s->input_used); return; }
        DEBUG_SERIAL_LINE_DECODER && fprintf(stderr,"sld:startbit_found at %d (%d)\n", startbit_start, iabs_samples_helper + startbit_start);

        //If the stop bit would be too far so that we reached the end of the buffer, then we return failed.
        //The caller can rearrange the buffer so that the whole character fits into it.
        float all_bits = 1 + s->databits + s->stopbits;
        DEBUG_SERIAL_LINE_DECODER && fprintf(stderr,"sld:all_bits = %f\n", all_bits);
        if(startbit_start + s->samples_per_bits * all_bits >= input_size) { s->input_used += MAX_M(0,startbit_start-2); DEBUG_SERIAL_LINE_DECODER && fprintf(stderr,"sld:return_stopbit_too_far (+%d)\n", s->input_used); return; }

        //We do the actual sampling.
        int di; //databit counter
        unsigned shr = 0;
        for(di=0; di < s->databits; di++)
        {
            int databit_start = startbit_start + (1+di+(0.5*(1-s->bit_sampling_width_ratio))) * s->samples_per_bits;
            int databit_end   = startbit_start + (1+di+(0.5*(1+s->bit_sampling_width_ratio))) * s->samples_per_bits;
            DEBUG_SERIAL_LINE_DECODER && fprintf(stderr,"sld:databit_start = %d (%d)\n", databit_start, iabs_samples_helper+databit_start);
            DEBUG_SERIAL_LINE_DECODER && fprintf(stderr,"sld:databit_end   = %d (%d)\n", databit_end,   iabs_samples_helper+databit_end);
            float databit_acc = 0;
            for(i=databit_start;i<databit_end;i++) { databit_acc += input[i]; /*DEBUG_SERIAL_LINE_DECODER && fprintf(stderr, "%f (%f) ", input[i], databit_acc);*/ }
            //DEBUG_SERIAL_LINE_DECODER && fprintf(stderr,"\n");
            DEBUG_SERIAL_LINE_DECODER && fprintf(stderr,"sld:databit_decision = %d\n", !!(databit_acc>0));
            shr=(shr<<1)|!!(databit_acc>0);
        }
        DEBUG_SERIAL_LINE_DECODER && fprintf(stderr,"sld:shr = 0x%x, %d\n", shr, shr);

        //We check if the stopbit is correct.
        int stopbit_start = startbit_start + (1+s->databits) * s->samples_per_bits + (s->stopbits * 0.5 * (1-s->bit_sampling_width_ratio)) * s->samples_per_bits;
        int stopbit_end   = startbit_start + (1+s->databits) * s->samples_per_bits + (s->stopbits * 0.5 * (1+s->bit_sampling_width_ratio)) * s->samples_per_bits;
        DEBUG_SERIAL_LINE_DECODER && fprintf(stderr,"sld:stopbit_start = %d (%d)\n", stopbit_start, iabs_samples_helper+stopbit_start);
        DEBUG_SERIAL_LINE_DECODER && fprintf(stderr,"sld:stopbit_end   = %d (%d)\n", stopbit_end,   iabs_samples_helper+stopbit_end);
        float stopbit_acc = 0;
        for(i=stopbit_start;i<stopbit_end;i++) { stopbit_acc += input[i]; DEBUG_SERIAL_LINE_DECODER && fprintf(stderr, "%f (%f) ", input[i], stopbit_acc); }
        DEBUG_SERIAL_LINE_DECODER && fprintf(stderr,"\n");
        if(stopbit_acc<0) { s->input_used += MIN_M(startbit_start + 1, input_size); DEBUG_SERIAL_LINE_DECODER && fprintf(stderr,"sld:return_stopbit_faulty (+%d)\n", s->input_used); return; }
        DEBUG_SERIAL_LINE_DECODER && fprintf(stderr,"sld:stopbit_found\n");

        //we write the output sample
        if(s->databits <= 8) output[s->output_size] = shr;
        else if(s->databits <= 16) output_s[s->output_size] = shr;
        else output_u[s->output_size] = shr;
        s->output_size++;

        int samples_used_up_now = MIN_M(startbit_start + all_bits * s->samples_per_bits, input_size);
        s->input_used += samples_used_up_now;
        input += samples_used_up_now;
        input_size -= samples_used_up_now;
        iabs_samples_helper += samples_used_up_now;
        if(!input_size) { DEBUG_SERIAL_LINE_DECODER && fprintf(stderr,"sld:return_no_more_input (+%d)\n", s->input_used); return; }
    }
    DEBUG_SERIAL_LINE_DECODER && fprintf(stderr, "sld: >> output_size = %d  (+%d)\n", s->output_size, s->input_used);
}

void generic_slicer_f_u8(float* input, unsigned char* output, int input_size, int n_symbols)
{
    float symbol_distance = 2.0/(n_symbols-1);
    for(int i=0;i<input_size;i++) 
        for(int j=0;j<n_symbols;j++)
        {
            float symbol_center = -1+j*symbol_distance;
            float symbol_low_limit = symbol_center-(symbol_distance/2);
            float symbol_high_limit = symbol_center+(symbol_distance/2);
            if(j==0)
            {
                if(input[i]<symbol_high_limit) 
                {
                    output[i]=j;
                    break;
                }
            }
            else if (j==n_symbols-1)
            {
                if(input[i]>=symbol_low_limit) 
                {
                    output[i]=j;
                    break;
                }
            }
            else 
            {
                if(input[i]>=symbol_low_limit && input[i]<symbol_high_limit) 
                {
                    output[i]=j;
                    break;
                }
            }
        }
}

void binary_slicer_f_u8(float* input, unsigned char* output, int input_size)
{
    for(int i=0;i<input_size;i++) output[i] = input[i] > 0;
}

void psk_modulator_u8_c(unsigned char* input, complexf* output, int input_size, int n_psk)
{
    //outputs one complex sample per input symbol
    float phase_increment = (2*M_PI)/n_psk;
    for(int i=0;i<input_size;i++)
    {
        float out_phase=phase_increment*input[i];
        iof(output,i)=cos(out_phase);
        qof(output,i)=sin(out_phase);
    }
}

void duplicate_samples_ntimes_u8_u8(unsigned char* input, unsigned char* output, int input_size_bytes, int sample_size_bytes, int ntimes)
{
    int l=0;
    for(int i=0;i<input_size_bytes;i+=sample_size_bytes)
        for(int k=0;k<ntimes;k++)
            for(int j=0;j<sample_size_bytes;j++)
                output[l++]=input[i+j];
}

complexf psk31_interpolate_sine_cc(complexf* input, complexf* output, int input_size, int interpolation, complexf last_input)
{
    int oi=0; //output index
    for(int i=0;i<input_size;i++)
    {
        for(int j=0; j<interpolation; j++)
        {
            float rate = (1+sin(-(M_PI/2)+M_PI*((j+1)/(float)interpolation)))/2;
            iof(output,oi)=iof(input,i) * rate + iof(&last_input,0) * (1-rate);
            qof(output,oi)=qof(input,i) * rate + qof(&last_input,0) * (1-rate);
            oi++;
        }
        last_input = input[i];
    }
    return last_input;
}

void pack_bits_1to8_u8_u8(unsigned char* input, unsigned char* output, int input_size)
{ //output size should be input_size × 8
    for(int i=0; i<input_size; i++)
        for(int bi=0; bi<8; bi++) //bi: bit index
            *(output++)=(input[i]>>bi)&1;
}


unsigned char pack_bits_8to1_u8_u8(unsigned char* input)
{
    unsigned char output;
    for(int i=0;i<8;i++)
    {
        output<<=1;
        output|=!!input[i];
    }
    return output;
}
unsigned char differential_codec(unsigned char* input, unsigned char* output, int input_size, int encode, unsigned char state)
{
    if(!encode)
        for(int i=0;i<input_size;i++) 
        {
            output[i] = input[i] == state;
            state = input[i];
        }
    else
        for(int i=0;i<input_size;i++) 
        {
            if(!input[i]) state=!state;
            output[i] = state;
        }
    return state;
}

/*
   _____                _                      _______ _           _                _____
  / ____|              (_)             ___    |__   __(_)         (_)              |  __ \
 | |     __ _ _ __ _ __ _  ___ _ __   ( _ )      | |   _ _ __ ___  _ _ __   __ _   | |__) |___  ___ _____   _____ _ __ _   _
 | |    / _` | '__| '__| |/ _ \ '__|  / _ \/\    | |  | | '_ ` _ \| | '_ \ / _` |  |  _  // _ \/ __/ _ \ \ / / _ \ '__| | | |
 | |___| (_| | |  | |  | |  __/ |    | (_>  <    | |  | | | | | | | | | | | (_| |  | | \ \  __/ (_| (_) \ V /  __/ |  | |_| |
  \_____\__,_|_|  |_|  |_|\___|_|     \___/\/    |_|  |_|_| |_| |_|_|_| |_|\__, |  |_|  \_\___|\___\___/ \_/ \___|_|   \__, |
                                                                            __/ |                                       __/ |
                                                                           |___/                                       |___/
*/

void pll_cc_init_pi_controller(pll_t* p, float bandwidth, float ko, float kd, float damping_factor)
{
    //kd: detector gain
    //ko: VCO gain
    float bandwidth_omega = 2*M_PI*bandwidth;
    p->alpha  = (damping_factor*2*bandwidth_omega)/(ko*kd);
    float sampling_rate = 1; //the bandwidth is normalized to the sampling rate
    p->beta   = (bandwidth_omega*bandwidth_omega)/(sampling_rate*ko*kd);
    p->iir_temp = p->dphase = p->output_phase = 0;
}

void pll_cc_init_p_controller(pll_t* p, float alpha)
{
    p->alpha = alpha;
    p->dphase=p->output_phase=0;
}


void pll_cc(pll_t* p, complexf* input, float* output_dphase, complexf* output_nco, int input_size)
{
    for(int i=0;i<input_size;i++)
    {
        p->output_phase += p->dphase;
        while(p->output_phase>PI) p->output_phase-=2*PI;
        while(p->output_phase<-PI) p->output_phase+=2*PI;
        complexf current_nco;
        iof(&current_nco,0) = sin(p->output_phase);
        qof(&current_nco,0) = cos(p->output_phase);
        if(output_nco) output_nco[i] = current_nco; //we don't output anything if it is a NULL pointer

        //accurate phase detector: calculating error from phase offset
        float input_phase = atan2(iof(input,i),qof(input,i));
        float new_dphase = input_phase - p->output_phase;
        while(new_dphase>PI) new_dphase-=2*PI;
        while(new_dphase<-PI) new_dphase+=2*PI;

        //modeling analog phase detector, which would be abs(input[i] * current_nco) if we had a real output signal, but what if we have complex signals?
        //qof(&current_nco,0)=-qof(&current_nco,0); //calculate conjugate
        //complexf multiply_result;
        //cmult(&multiply_result, &input[i], &current_nco);
        //output_nco[i] = multiply_result;
        //float new_dphase = absof(&multiply_result,0);

        if(p->pll_type == PLL_PI_CONTROLLER)
        {
            p->dphase = new_dphase * p->alpha + p->iir_temp;
            p->iir_temp += new_dphase * p->beta;

            while(p->dphase>PI) p->dphase-=2*PI; //won't need this one
            while(p->dphase<-PI) p->dphase+=2*PI;
        }
        else if(p->pll_type == PLL_P_CONTROLLER)
        {
            p->dphase = new_dphase * p->alpha;
        }
        else return;
        if(output_dphase) output_dphase[i] = -p->dphase;
        //if(output_dphase) output_dphase[i] = new_dphase/10;
    }
}

void octave_plot_point_on_cplxsig(complexf* signal, int signal_size, float error, int index, int correction_offset, char* writefiles_path, int points_size, ...)
{
    static int figure_output_counter = 0;
    int* points_z = (int*)malloc(sizeof(int)*points_size);
    int* points_color = (int*)malloc(sizeof(int)*points_size);
    va_list vl;
    va_start(vl,points_size);
    for(int i=0;i<points_size;i++)
    {
        points_z[i] = va_arg(vl, int);
        points_color[i] = va_arg(vl, int);
    }
    if(writefiles_path && !figure_output_counter) fprintf(stderr, "cf=figure();\n");
    fprintf(stderr, "N = %d;\nisig = [", signal_size);
    for(int i=0;i<signal_size;i++) fprintf(stderr, "%f ", iof(signal,i));
    fprintf(stderr, "];\nqsig = [");
    for(int i=0;i<signal_size;i++) fprintf(stderr, "%f ", qof(signal,i));
    fprintf(stderr, "];\nzsig = [0:N-1];\nsubplot(2,2,[2 4]);\nplot3(isig,zsig,qsig,\"b-\",");
    for(int i=0;i<points_size;i++)
        fprintf(stderr, "[%f],[%d],[%f],\"%c.\"%c", 
            iof(signal, points_z[i]), points_z[i], qof(signal, points_z[i]), 
            (char)points_color[i]&0xff, (i<points_size-1)?',':' '
        );
    va_end(vl);
    fprintf(stderr, ");\ntitle(\"index = %d, error = %f, cxoffs = %d\");\nsubplot(2,2,1);\nplot(zsig, isig,\"b-\",", index, error, correction_offset);
    for(int i=0;i<points_size;i++)
        fprintf(stderr, "[%d],[%f],\"%c.\"%c", 
            points_z[i], iof(signal, points_z[i]),
            (char)points_color[i]&0xff, (i<points_size-1)?',':' '
        );
    fprintf(stderr, ");\nsubplot(2,2,3);\nplot(zsig, qsig,\"b-\",");
    for(int i=0;i<points_size;i++)
        fprintf(stderr, "[%d],[%f],\"%c.\"%c", 
            points_z[i], qof(signal, points_z[i]),
            (char)points_color[i]&0xff, (i<points_size-1)?',':' '
        );
    fprintf(stderr, ");\n");
    if(writefiles_path) fprintf(stderr, "print(cf, \"%s/%05d.png\", \"-S1024,1024\");\n", writefiles_path, figure_output_counter++); 
    fflush(stderr);
    free(points_z);
    free(points_color);
}

timing_recovery_state_t timing_recovery_init(timing_recovery_algorithm_t algorithm, int decimation_rate, int use_q, float loop_gain, float max_error, int debug_every_nth, char* debug_writefiles_path)
{
    timing_recovery_state_t to_return;
    to_return.algorithm = algorithm;
    to_return.decimation_rate = decimation_rate;
    to_return.loop_gain = loop_gain;
    to_return.max_error = max_error;
    to_return.use_q = use_q;
    to_return.debug_phase = to_return.debug_every_nth = debug_every_nth; //debug is effective if it is >=0
    to_return.last_correction_offset = 0;
    to_return.earlylate_ratio = 0.25; //0..0.5
    to_return.debug_writefiles_path = debug_writefiles_path;
    return to_return;
}

#define MTIMINGR_HDEBUG 0

void timing_recovery_cc(complexf* input, complexf* output, int input_size, float* timing_error, int* sampled_indexes, timing_recovery_state_t* state)
{
    //We always assume that the input starts at center of the first symbol cross before the first symbol.
    //Last time we consumed that much from the input samples that it is there.
    int correction_offset = state->last_correction_offset;
    int current_bitstart_index = 0;
    int num_samples_bit = state->decimation_rate;
    int num_samples_halfbit = state->decimation_rate / 2;
    int num_samples_quarterbit = state->decimation_rate / 4;
    int num_samples_earlylate_wing = num_samples_bit * state->earlylate_ratio;
    float error;
    int el_point_left_index, el_point_right_index, el_point_mid_index;
    int si = 0;
    if(state->debug_every_nth>=0) fprintf(stderr, "disp(\"begin timing_recovery_cc\");\n");
    if(MTIMINGR_HDEBUG) fprintf(stderr, "timing_recovery_cc started, nsb = %d, nshb = %d, nsqb = %d\n", num_samples_bit, num_samples_halfbit, num_samples_quarterbit);
    {
        for(;;)
        {
            //the MathWorks style algorithm has correction_offset.
            //correction_offset = 0;            
            if(current_bitstart_index + num_samples_halfbit * 3 >= input_size) break;
            if(MTIMINGR_HDEBUG) fprintf(stderr, "current_bitstart_index = %d, input_size = %d, correction_offset(prev) = %d\n", 
                    current_bitstart_index, input_size, correction_offset);
            
            if(correction_offset<=-num_samples_quarterbit*0.9 || correction_offset>=0.9*num_samples_quarterbit)
            {
                if(MTIMINGR_HDEBUG) fprintf(stderr, "correction_offset = %d, reset to 0!\n", correction_offset); 
                correction_offset = 0;
            }
            //should check if the sign of the correction_offset (or disabling it) has an effect on the EVM.
            //it is also a possibility to disable multiplying with the magnitude
            if(state->algorithm == TIMING_RECOVERY_ALGORITHM_EARLYLATE)
            {
                //bitstart index should be at symbol edge, maximum effect point is at current_bitstart_index + num_samples_halfbit
                el_point_right_index  = current_bitstart_index + num_samples_earlylate_wing * 3;
                el_point_left_index   = current_bitstart_index + num_samples_earlylate_wing * 1 - correction_offset;
                el_point_mid_index    = current_bitstart_index + num_samples_halfbit;
                if(sampled_indexes) sampled_indexes[si]=el_point_mid_index;
                output[si++] = input[el_point_mid_index];
            }
            else if(state->algorithm == TIMING_RECOVERY_ALGORITHM_GARDNER)
            {
                //maximum effect point is at current_bitstart_index
                el_point_right_index  = current_bitstart_index + num_samples_halfbit * 3;
                el_point_left_index   = current_bitstart_index + num_samples_halfbit * 1;
                el_point_mid_index    = current_bitstart_index + num_samples_halfbit * 2;
                if(sampled_indexes) sampled_indexes[si]=el_point_left_index;
                output[si++] = input[el_point_left_index];
            }
            else break;

            error = ( iof(input, el_point_right_index) - iof(input, el_point_left_index) ) * iof(input, el_point_mid_index); 
            if(state->use_q)
            {
                error += ( qof(input, el_point_right_index) - qof(input, el_point_left_index)) * qof(input, el_point_mid_index); 
                error /= 2;
            }
            //Original correction method: this version can only move a single sample in any direction
            //current_bitstart_index += num_samples_halfbit * 2 + (error)?((error<0)?1:-1):0;

            if(timing_error) timing_error[si-1]=error; //it is not written if NULL
            
            if(error>state->max_error) error=state->max_error;
            if(error<-state->max_error) error=-state->max_error;
            if(state->debug_every_nth>=0)
            {
                if(state->debug_every_nth==0 || state->debug_phase==0) 
                {
                    state->debug_phase = state->debug_every_nth;
                    octave_plot_point_on_cplxsig(input+current_bitstart_index, state->decimation_rate*2, 
                        error, 
                        current_bitstart_index, 
                        correction_offset,
                        state->debug_writefiles_path,
                        3,
                        el_point_left_index - current_bitstart_index,  'r',
                        el_point_right_index - current_bitstart_index, 'r',
                        el_point_mid_index - current_bitstart_index,   'r',
                        0);
                }
                else state->debug_phase--;
            }
            int error_sign = (state->algorithm == TIMING_RECOVERY_ALGORITHM_GARDNER) ? -1 : 1;
            correction_offset = num_samples_halfbit * error_sign * error * state->loop_gain;
            current_bitstart_index += num_samples_bit + correction_offset;
            if(si>=input_size) 
            { 
                if(MTIMINGR_HDEBUG) fprintf(stderr, "oops_out_of_si!\n"); 
                break; 
            }
        }
    }
    state->input_processed = current_bitstart_index;
    state->output_size = si;
    state->last_correction_offset = correction_offset;
}

#define MTIMINGR_GAS(NAME) \
    if(!strcmp( #NAME , input )) return TIMING_RECOVERY_ALGORITHM_ ## NAME;

timing_recovery_algorithm_t timing_recovery_get_algorithm_from_string(char* input)
{
    MTIMINGR_GAS(GARDNER);
    MTIMINGR_GAS(EARLYLATE);
    return TIMING_RECOVERY_ALGORITHM_DEFAULT;
}

#define MTIMINGR_GSA(NAME) \
    if(algorithm == TIMING_RECOVERY_ALGORITHM_ ## NAME) return #NAME;

char* timing_recovery_get_string_from_algorithm(timing_recovery_algorithm_t algorithm)
{
    MTIMINGR_GSA(GARDNER);
    MTIMINGR_GSA(EARLYLATE);
    return "INVALID";
}

void init_bpsk_costas_loop_cc(bpsk_costas_loop_state_t* s, int decision_directed, float damping_factor, float bandwidth)
{
    //fprintf(stderr, "init_bpsk_costas_loop_cc: bandwidth = %f, damping_factor = %f\n", bandwidth, damping_factor);
    //based on: http://gnuradio.squarespace.com/blog/2011/8/13/control-loop-gain-values.html
    float bandwidth_omega = 2*PI*bandwidth; //so that the bandwidth should be around 0.01 by default (2pi/100), and the damping_factor should be default 0.707
    float denomiator = 1+2*damping_factor*bandwidth_omega+bandwidth_omega*bandwidth_omega;
    fprintf(stderr, "damp = %f, bw = %f, bwomega = %f\n", damping_factor, bandwidth, bandwidth_omega);
    s->alpha = (4*damping_factor*bandwidth_omega)/denomiator;
    s->beta = (4*bandwidth_omega*bandwidth_omega)/denomiator;
    s->current_freq = s->dphase = s->nco_phase = 0;
    s->dphase_max=bandwidth_omega; //this has been determined by experiment: if dphase is out of [-dphase_max; dphase_max] it might actually hang and not come back 
    s->dphase_max_reset_to_zero=0;
}

void bpsk_costas_loop_cc(complexf* input, complexf* output, int input_size, float* output_error, float* output_dphase, complexf* output_nco, bpsk_costas_loop_state_t* s)
{
    for(int i=0;i<input_size;i++)
    {
        complexf nco_sample;
        e_powj(&nco_sample, s->nco_phase);
        cmult(&output[i], &input[i], &nco_sample);
        if(output_nco) output_nco[i]=nco_sample;
        float error = 0;
        if(s->decision_directed)
        {
            float output_phase = atan2(qof(output,i),iof(output,i));
            if (fabs(output_phase)<PI/2) 
                error = -output_phase;
            else
            {
                error = PI-output_phase;
                while(error>PI) error -= 2*PI;
            }
        }
        else error = PI*iof(output,i)*qof(output,i);
        if(output_error) output_error[i]=error;
        s->current_freq += error * s->beta;
        s->dphase = error * s->alpha + s->current_freq;
        if(s->dphase>s->dphase_max)  s->dphase = (s->dphase_max_reset_to_zero) ? 0 : s->dphase_max;
        if(s->dphase<-s->dphase_max) s->dphase = (s->dphase_max_reset_to_zero) ? 0 : -s->dphase_max;
        if(output_dphase) output_dphase[i]=s->dphase;
        //fprintf(stderr, "  error = %f; dphase = %f; nco_phase = %f;\n", error, s->dphase, s->nco_phase);

        //step NCO
        s->nco_phase += s->dphase;
        while(s->nco_phase>2*PI) s->nco_phase-=2*PI;
        while(s->nco_phase<=0) s->nco_phase+=2*PI;
    }
}

#if 0
bpsk_costas_loop_state_t init_bpsk_costas_loop_cc(float samples_per_bits) 
{ 
    bpsk_costas_loop_state_t state;
    state.vco_phase = 0;
    state.last_vco_phase_addition = 0;
    float virtual_sampling_rate = 10000;
    float virtual_data_rate = virtual_sampling_rate / samples_per_bits;
    fprintf(stderr, "virtual_sampling_rate = %g, virtual_data_rate = %g\n", virtual_sampling_rate, virtual_data_rate);
    float rc_filter_cutoff = virtual_data_rate * 2; //this is so far the best
    float rc_filter_rc = 1/(2*M_PI*rc_filter_cutoff); //as of Equation 24 in Feigin
    float virtual_sampling_dt = 1.0/virtual_sampling_rate;
    fprintf(stderr, "rc_filter_cutoff = %g, rc_filter_rc = %g, virtual_sampling_dt = %g\n", 
        rc_filter_cutoff, rc_filter_rc, virtual_sampling_dt);
    state.rc_filter_alpha = virtual_sampling_dt/(rc_filter_rc+virtual_sampling_dt); //https://en.wikipedia.org/wiki/Low-pass_filter
    float rc_filter_omega_cutoff = 2*M_PI*rc_filter_cutoff;
    state.vco_phase_addition_multiplier = 8*rc_filter_omega_cutoff / (virtual_sampling_rate); //as of Equation 25 in Feigin, assuming input signal amplitude of 1 (to 1V) and (state.vco_phase_addition_multiplier*<vco_input>), a value in radians, will be added to the vco_phase directly.
    fprintf(stderr, "rc_filter_alpha = %g, rc_filter_omega_cutoff = %g, vco_phase_addition_multiplier = %g\n", 
            state.rc_filter_alpha, rc_filter_omega_cutoff, state.vco_phase_addition_multiplier);
    return state;
}

void bpsk_costas_loop_c1mc(complexf* input, complexf* output, int input_size, bpsk_costas_loop_state_t* state)
{
    int debug = 0;
    if(debug) fprintf(stderr, "costas:\n");
    for(int i=0;i<input_size;i++)
    {
        float input_phase = atan2(input[i].q, input[i].i);
        float input_and_vco_mixed_phase = input_phase - state->vco_phase;
        if(debug) fprintf(stderr, "%g | %g\n", input_and_vco_mixed_phase, input_phase), debug--;
        complexf input_and_vco_mixed_sample; 
        e_powj(&input_and_vco_mixed_sample, input_and_vco_mixed_phase);
        
        complexf vco_sample;
        e_powj(&vco_sample, -state->vco_phase);
        //cmult(&input_and_vco_mixed_sample, &input[i], &vco_sample);//if this is enabled, the real input sample is used, not the amplitude normalized 

        float loop_output_i = 
            input_and_vco_mixed_sample.i * state->rc_filter_alpha + state->last_lpfi_output * (1-state->rc_filter_alpha);
        float loop_output_q = 
            input_and_vco_mixed_sample.q * state->rc_filter_alpha + state->last_lpfq_output * (1-state->rc_filter_alpha);
        //loop_output_i = input_and_vco_mixed_sample.i;
        //loop_output_q = input_and_vco_mixed_sample.q;
        state->last_lpfi_output = loop_output_i;
        state->last_lpfq_output = loop_output_q;
        float vco_phase_addition = loop_output_i * loop_output_q * state->vco_phase_addition_multiplier;
        //vco_phase_addition = vco_phase_addition * state->rc_filter_alpha + state->last_vco_phase_addition * (1-state->rc_filter_alpha);
        //state->last_vco_phase_addition = vco_phase_addition;
        state->vco_phase += vco_phase_addition;
        while(state->vco_phase>PI) state->vco_phase-=2*PI;
        while(state->vco_phase<-PI) state->vco_phase+=2*PI;
        cmult(&output[i], &input[i], &vco_sample);
    }
}
#endif

void simple_agc_cc(complexf* input, complexf* output, int input_size, float rate, float reference, float max_gain, float* current_gain)
{
    float rate_1minus=1-rate;
    int debugn = 0;
    for(int i=0;i<input_size;i++)
    {
        float amplitude = sqrt(input[i].i*input[i].i+input[i].q*input[i].q);
        float ideal_gain = (reference/amplitude);
        if(ideal_gain>max_gain) ideal_gain = max_gain;
        if(ideal_gain<=0) ideal_gain = 0;
        //*current_gain += (ideal_gain-(*current_gain))*rate;
        *current_gain = (ideal_gain-(*current_gain))*rate + (*current_gain)*rate_1minus;
        //if(debugn<100) fprintf(stderr, "cgain: %g\n", *current_gain), debugn++;
        output[i].i=(*current_gain)*input[i].i;
        output[i].q=(*current_gain)*input[i].q;
    }
}

void firdes_add_peak_c(complexf* output, int length, float rate, window_t window, int add, int normalize)
{
    //add=0: malloc output previously
    //add=1: calloc output previously
    complexf* taps = (complexf*)malloc(sizeof(complexf)*length);
    int middle=length/2;
    float phase = 0, phase_addition = -rate*M_PI*2;
    float (*window_function)(float) = firdes_get_window_kernel(window);
    for(int i=0; i<length; i++) //@@firdes_add_peak_c: calculate taps
    {
        e_powj(&taps[i], phase);
        float window_multiplier = window_function(fabs((float)(middle-i)/middle));
        taps[i].i *= window_multiplier;
        taps[i].q *= window_multiplier;
        phase += phase_addition;
        while(phase>2*M_PI) phase-=2*M_PI;
        while(phase<0) phase+=2*M_PI;
    }

    //Normalize filter kernel
    if(add) 
        for(int i=0;i<length;i++)
        {
            output[i].i += taps[i].i;
            output[i].q += taps[i].q;
        }
    else for(int i=0;i<length;i++) output[i] = taps[i];
    if(normalize)
    {
        float sum=0;
        for(int i=0;i<length;i++) //@firdes_add_peak_c: normalize pass 1
        {
            sum+=sqrt(output[i].i*output[i].i + output[i].q*output[i].q);
        }
        for(int i=0;i<length;i++) //@firdes_add_peak_c: normalize pass 2
        {
            output[i].i/=sum;
            output[i].q/=sum;
        }
    }
}

int apply_fir_cc(complexf* input, complexf* output, int input_size, complexf* taps, int taps_length)
{
    int i;
    for(i=0; i<input_size-taps_length+1; i++)
    {
        csetnull(&output[i]);
        for(int ti=0;ti<taps_length;ti++)
        {
            cmultadd(&output[i], &input[i+ti], &taps[ti]);
        }
    }
    return i;
}


int apply_real_fir_cc(complexf* input, complexf* output, int input_size, float* taps, int taps_length)
{
    int i;
    for(i=0; i<input_size-taps_length+1; i++)
    {
        float acci = 0, accq = 0;
        for(int ti=0;ti<taps_length;ti++)
        {
            acci += iof(input,i+ti)*taps[ti];
            accq += qof(input,i+ti)*taps[ti];
        }
        iof(output,i)=acci;
        qof(output,i)=accq;
    }
    return i;
}

float normalized_timing_variance_u32_f(unsigned* input, float* temp, int input_size, int samples_per_symbol, int initial_sample_offset, int debug_print)
{
    float *ndiff_rad = temp;
    float ndiff_rad_mean = 0;
    for(int i=0;i<input_size;i++) 
    {
        //find out which real sample index this input sample index is the nearest to.
        unsigned sinearest = (input[i]-initial_sample_offset) / samples_per_symbol;
        unsigned sinearest_remain = (input[i]-initial_sample_offset) % samples_per_symbol;
        if(sinearest_remain>samples_per_symbol/2) sinearest++;
        unsigned socorrect = initial_sample_offset+(sinearest*samples_per_symbol); //the sample offset which input[i] should have been, in order to sample at the maximum effect point
        int sodiff = abs(socorrect-input[i]);
        float ndiff = (float)sodiff/samples_per_symbol;

        ndiff_rad[i] = ndiff*PI;
        ndiff_rad_mean = ndiff_rad_mean*(((float)i)/(i+1))+(ndiff_rad[i]/(i+1));
        if(debug_print) fprintf(stderr, "input[%d] = %u, sinearest = %u, socorrect = %u, sodiff = %u, ndiff = %f, ndiff_rad[i] = %f, ndiff_rad_mean = %f\n", i, input[i], sinearest, socorrect, sodiff, ndiff, ndiff_rad[i], ndiff_rad_mean);
    }
    fprintf(stderr, "ndiff_rad_mean = %f\n", ndiff_rad_mean);

    float result = 0;
    for(int i=0;i<input_size;i++) result+=(powf(ndiff_rad[i]-ndiff_rad_mean,2))/(input_size-1);
    //fprintf(stderr, "nv = %f\n", result);
    return result;
}

void dbpsk_decoder_c_u8(complexf* input, unsigned char* output, int input_size)
{
    static complexf last_input;
    for(int i=0;i<input_size;i++)
    {
        float phase = atan2(qof(input,i), iof(input,i));
        float last_phase = atan2(qofv(last_input), iofv(last_input));
        float dphase = phase-last_phase;
        while(dphase<-PI) dphase+=2*PI;
        while(dphase>=PI) dphase-=2*PI;
        if( (dphase>(PI/2)) || (dphase<(-PI/2)) ) output[i]=0;
        else output[i]=1;
        last_input = input[i];
    }
}

int bfsk_demod_cf(complexf* input, float* output, int input_size, complexf* mark_filter, complexf* space_filter, int taps_length)
{
    complexf acc_space, acc_mark;
    for(int i=0; i<input_size-taps_length+1; i++)
    {
        csetnull(&acc_space);
        csetnull(&acc_mark);
        for(int ti=0;ti<taps_length;ti++)
        {
            cmultadd(&acc_space, &input[i+ti], &space_filter[ti]);
            cmultadd(&acc_mark,  &input[i+ti], &mark_filter[ti]);
        }
        output[i] = - ( iofv(acc_space)*iofv(acc_space)+qofv(acc_space)*qofv(acc_space) ) +
            ( iofv(acc_mark)*iofv(acc_mark)+qofv(acc_mark)*qofv(acc_mark) );
    }
    return input_size-taps_length+1;
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

void convert_s8_f(signed char* input, float* output, int input_size)
{
    for(int i=0;i<input_size;i++) output[i]=((float)input[i])/SCHAR_MAX; //@convert_s8_f
}

void convert_s16_f(short* input, float* output, int input_size)
{
    for(int i=0;i<input_size;i++) output[i]=(float)input[i]/SHRT_MAX; //@convert_s16_f
}

void convert_f_u8(float* input, unsigned char* output, int input_size)
{
    for(int i=0;i<input_size;i++) output[i]=input[i]*UCHAR_MAX*0.5+128; //@convert_f_u8
    //128 above is the correct value to add. In any other case a DC component
    //of at least -60 dB is shown on the FFT plot after convert_f_u8 -> convert_u8_f
}

void convert_f_s8(float* input, signed char* output, int input_size)
{
    for(int i=0;i<input_size;i++) output[i]=input[i]*SCHAR_MAX; //@convert_f_s8
}

void convert_f_s16(float* input, short* output, int input_size)
{
    /*for(int i=0;i<input_size;i++)
    {
        if(input[i]>1.0) input[i]=1.0;
        if(input[i]<-1.0) input[i]=-1.0;
    }*/
    for(int i=0;i<input_size;i++) output[i]=input[i]*SHRT_MAX; //@convert_f_s16
}

void convert_i16_f(short* input, float* output, int input_size) { convert_s16_f(input, output, input_size); }
void convert_f_i16(float* input, short* output, int input_size) { convert_f_s16(input, output, input_size); }

void convert_f_s24(float* input, unsigned char* output, int input_size, int bigendian)
{
    int k=0;
    if(bigendian) for(int i=0;i<input_size;i++)
    {
        int temp=input[i]*(INT_MAX>>8);
        unsigned char* ptemp=(unsigned char*)&temp;
        output[k++]=*ptemp;
        output[k++]=*(ptemp+1);
        output[k++]=*(ptemp+2);
    }
    else for(int i=0;i<input_size;i++)
    {
        int temp=input[i]*(INT_MAX>>8);
        unsigned char* ptemp=(unsigned char*)&temp;
        output[k++]=*(ptemp+2);
        output[k++]=*(ptemp+1);
        output[k++]=*ptemp;
    }
}

void convert_s24_f(unsigned char* input, float* output, int input_size, int bigendian)
{
    int k=0;
    if(bigendian) for(int i=0;i<input_size*3;i+=3)
    {
        int temp=(input[i+2]<<24)|(input[i+1]<<16)|(input[i]<<8);
        output[k++]=temp/(float)(INT_MAX-256);
    }
    else for(int i=0;i<input_size*3;i+=3)
    {
        int temp=(input[i+2]<<8)|(input[i+1]<<16)|(input[i]<<24);
        output[k++]=temp/(float)(INT_MAX-256);
    }
}

FILE* init_get_random_samples_f()
{
    return fopen("/dev/urandom", "r");
}

void get_random_samples_f(float* output, int output_size, FILE* status)
{
    int* pioutput = (int*)output;
    fread((unsigned char*)output, sizeof(float), output_size, status);
    for(int i=0;i<output_size;i++)
    {
        float tempi = pioutput[i];
        output[i] = tempi/((float)(INT_MAX)); //*0.82
    }
}

void get_random_gaussian_samples_c(complexf* output, int output_size, FILE* status)
{
    int* pioutput = (int*)output;
    fread((unsigned char*)output, sizeof(complexf), output_size, status);
    for(int i=0;i<output_size;i++)
    {
        float u1 = 0.5+0.49999999*(((float)pioutput[2*i])/(float)INT_MAX);
        float u2 = 0.5+0.49999999*(((float)pioutput[2*i+1])/(float)INT_MAX);
        iof(output, i)=sqrt(-2*log(u1))*cos(2*PI*u2);
        qof(output, i)=sqrt(-2*log(u1))*sin(2*PI*u2);
    }
}

int deinit_get_random_samples_f(FILE* status)
{
    return fclose(status);
}

int firdes_cosine_f(float* taps, int taps_length, int samples_per_symbol)
{
    //needs a taps_length 2 × samples_per_symbol + 1
    int middle_i=taps_length/2;
    for(int i=0;i<samples_per_symbol;i++) taps[middle_i+i]=taps[middle_i-i]=(1+cos(PI*i/(float)samples_per_symbol))/2;
    //for(int i=0;i<taps_length;i++) taps[i]=powf(taps[i],2);
    normalize_fir_f(taps, taps, taps_length);
}

int firdes_rrc_f(float* taps, int taps_length, int samples_per_symbol, float beta)
{
    //needs an odd taps_length
    int middle_i=taps_length/2;
    taps[middle_i]=(1/(float)samples_per_symbol)*(1+beta*(4/PI-1));
    for(int i=1;i<1+taps_length/2;i++) 
    {
        if(i==samples_per_symbol/(4*beta)) 
            taps[middle_i+i]=taps[middle_i-i]=(beta/(samples_per_symbol*sqrt(2)))*((1+(2/PI))*sin(PI/(4*beta))+(1-(2/PI))*cos(PI/(4*beta)));
        else
            taps[middle_i+i]=taps[middle_i-i]=(1/(float)samples_per_symbol)*
                (sin(PI*(i/(float)samples_per_symbol)*(1-beta)) + 4*beta*(i/(float)samples_per_symbol)*cos(PI*(i/(float)samples_per_symbol)*(1+beta)))/
                (PI*(i/(float)samples_per_symbol)*(1-powf(4*beta*(i/(float)samples_per_symbol),2)));
    }
    normalize_fir_f(taps, taps, taps_length);
}

void plain_interpolate_cc(complexf* input, complexf* output, int input_size, int interpolation)
{
    for(int i=0;i<input_size;i++)
    {
        output[i*interpolation]=input[i];
        bzero(output+(interpolation*i)+1, (interpolation-1)*sizeof(complexf));
    }
}

#define MMATCHEDFILT_GAS(NAME) \
    if(!strcmp( #NAME , input )) return MATCHED_FILTER_ ## NAME;

matched_filter_type_t matched_filter_get_type_from_string(char* input)
{
    MMATCHEDFILT_GAS(RRC);
    MMATCHEDFILT_GAS(COSINE);
    return MATCHED_FILTER_DEFAULT;
}

float* add_ff(float* input1, float* input2, float* output, int input_size)
{
    for(int i=0;i<input_size;i++) output[i]=input1[i]+input2[i];
}

float* add_const_cc(complexf* input, complexf* output, int input_size, complexf x)
{
    for(int i=0;i<input_size;i++)
    {
        iof(output,i)=iof(input,i)+iofv(x);
        qof(output,i)=iof(input,i)+qofv(x);
    }
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
