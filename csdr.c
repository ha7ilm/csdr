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

#define _POSIX_C_SOURCE 199309L
#define _BSD_SOURCE
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <time.h>
#include <stdarg.h>
#include "libcsdr.h"
#include "libcsdr_gpl.h"
#include "ima_adpcm.h"
#include <sched.h>
#include <math.h>

char usage[]=
"csdr - a simple commandline tool for Software Defined Radio receiver DSP.\n\n"
"usage: \n\n"
"    csdr function_name <function_param1> <function_param2> [optional_param] ...\n\n"
"list of functions:\n\n"
"    convert_u8_f\n"
"    convert_f_u8\n"
"    convert_s8_f\n"
"    convert_f_s8\n"
"    convert_f_s16\n"
"    convert_s16_f\n"
"    realpart_cf\n"
"    clipdetect_ff\n"
"    limit_ff [max_amplitude]\n"
"    gain_ff <gain>\n"
"    clone\n"
"    none\n"
"    yes_f <to_repeat> [buf_times]\n"
"    detect_nan_ff\n"
"    floatdump_f\n"
"    flowcontrol <data_rate> <reads_per_second> [prebuffer_sec] [thrust]\n"
"    shift_math_cc <rate>\n"
"    shift_math_cc --fifo <fifo_path>\n"
"    shift_addition_cc <rate>\n"
"    shift_addition_cc --fifo <fifo_path>\n"
"    shift_addition_cc_test\n"
"    shift_table_cc <rate> [table_size]\n"
"    decimating_shift_addition_cc <rate> [decimation]\n"
"    dcblock_ff\n"
"    fastdcblock_ff\n"
"    fmdemod_atan_cf\n"
"    fmdemod_quadri_cf\n"
"    fmdemod_quadri_novect_cf\n"
"    deemphasis_wfm_ff <sample_rate> <tau>\n"
"    deemphasis_nfm_ff <one_of_the_predefined_sample_rates>\n"
"    amdemod_cf\n"
"    amdemod_estimator_cf\n"
"    fir_decimate_cc <decimation_factor> [transition_bw [window]]\n"
"    firdes_lowpass_f <cutoff_rate> <length> [window [--octave]]\n"
"    firdes_bandpass_c <low_cut> <high_cut> <length> [window [--octave]]\n"
"    agc_ff [hang_time [reference [attack_rate [decay_rate [max_gain [attack_wait [filter_alpha]]]]]]]\n"
"    fastagc_ff [block_size [reference]]\n"
"    rational_resampler_ff <interpolation> <decimation> [transition_bw [window]]\n"
"    fractional_decimator_ff <decimation_rate> [transition_bw [window]]\n"
"    fft_cc <fft_size> <out_of_every_n_samples> [window [--octave] [--benchmark]]\n"
"    logpower_cf [add_db]\n"
"    fft_benchmark <fft_size> <fft_cycles> [--benchmark]\n"
"    bandpass_fir_fft_cc <low_cut> <high_cut> <transition_bw> [window]\n"
"    bandpass_fir_fft_cc --fifo <fifo_path> <transition_bw> [window]\n"
"    encode_ima_adpcm_s16_u8\n"
"    decode_ima_adpcm_u8_s16\n"
"    compress_fft_adpcm_f_u8 <fft_size>\n"
"    flowcontrol <data_rate> <reads_per_second>\n"
"    through\n"
"    dsb_fc [q_value]\n"
"    convert_f_samperf <wait_for_this_sample> \n"
"    fmmod_fc\n"
"    fixed_amplitude_cc <new_amplitude>\n"
"    monos2stereo_s16\n"
"    setbuf <buffer_size>\n"
"    fft_exchange_sides_ff <fft_size>\n"
"    squelch_and_smeter_cc --fifo <squelch_fifo> --outfifo <smeter_fifo> <use_every_nth> <report_every_nth>\n"
"    \n"
;

//change on 2015-08-29: we rather dynamically determine the bufsize
//#define BUFSIZE (1024)
//#define BIG_BUFSIZE (1024*16)
//should be multiple of 16! (size of double complex)
//also, keep in mind that shift_addition_cc works better the smaller this buffer is.

int env_csdr_fixed_bufsize = 1024;
int env_csdr_fixed_big_bufsize = 1024*16;
int env_csdr_dynamic_bufsize_on = 0;
int env_csdr_print_bufsizes = 0;
int bigbufs = 0;

//change on on 2015-08-29: we don't yield at all. fread() will do it if it blocks
#define YIELD_EVERY_N_TIMES 3
//#define TRY_YIELD if(++yield_counter%YIELD_EVERY_N_TIMES==0) sched_yield()
#define TRY_YIELD
//unsigned yield_counter=0;

int badsyntax(char* why)
{
	if(why==0) fprintf(stderr, "%s", usage);
	else fprintf(stderr, "csdr: %s\n\n", why);
	return -1;
}

int clipdetect_ff(float* input, int input_size)
{
	for(int i=0;i<input_size;i++)
	{
		if(input[i]<-1.0) { fprintf(stderr, "csdr clipdetect_ff: Signal value below -1.0!\n"); return -1; }
		if(input[i]>1.0) { fprintf(stderr, "csdr clipdetect_ff: Signal value above 1.0!\n"); return 1; }
	}
	return 0;
}

int clone_(int bufsize_param)
{
		unsigned char* clone_buffer;
		clone_buffer = (unsigned char*)malloc(bufsize_param*sizeof(unsigned char));
		for(;;)
		{
			fread(clone_buffer, sizeof(unsigned char), bufsize_param, stdin);
			fwrite(clone_buffer, sizeof(unsigned char), bufsize_param, stdout);
			TRY_YIELD;
		}
}

#define FREAD_R		fread (input_buffer,	sizeof(float),		the_bufsize, stdin)
#define FREAD_C		fread (input_buffer,	sizeof(float)*2,	the_bufsize, stdin)
#define FWRITE_R 	fwrite (output_buffer,	sizeof(float),		the_bufsize, stdout)
#define FWRITE_C	fwrite (output_buffer,	sizeof(float)*2,	the_bufsize, stdout)
#define FEOF_CHECK	if(feof(stdin)) return 0
//#define BIG_FREAD_C fread(input_buffer, sizeof(float)*2, BIG_BUFSIZE, stdin)
//#define BIG_FWRITE_C fwrite(output_buffer, sizeof(float)*2, BIG_BUFSIZE, stdout)

int init_fifo(int argc, char *argv[])
{
	if(argc>=4)
	{
		if(!strcmp(argv[2],"--fifo"))
		{
			fprintf(stderr,"csdr: fifo control mode on\n");
			int fd = open(argv[3], O_RDONLY);
			int flags = fcntl(fd, F_GETFL, 0);
			fcntl(fd, F_SETFL, flags | O_NONBLOCK);
			return fd;
		}
	}
	return 0;
}



#define RFCTL_BUFSIZE 1024

int read_fifo_ctl(int fd, char* format, ...)
{
	if(!fd) return 0;
	static char buffer[RFCTL_BUFSIZE];
	static int buffer_index=0;
	int bytes_read=read(fd,buffer+buffer_index,(RFCTL_BUFSIZE-buffer_index)*sizeof(char));
	if(bytes_read<=0) return 0;

	int prev_newline_at=0;
	int last_newline_at=0;
	for(int i=0;i<buffer_index+bytes_read;i++)
	{
		if(buffer[i]=='\n')
		{
			prev_newline_at=last_newline_at;
			last_newline_at=i+1;
		}
	}
	if(last_newline_at)
	{
		//fprintf(stderr,"pna=%d lna=%d\n",prev_newline_at,last_newline_at);
		va_list vl;
		va_start(vl,format);
		vsscanf(buffer+prev_newline_at,format,vl);
		va_end(vl);
		memmove(buffer,buffer+last_newline_at,buffer_index+bytes_read-last_newline_at);
		buffer_index=bytes_read-last_newline_at;
		return 1;
	}
	else
	{
		buffer_index+=bytes_read;
	 	return 0;
	}
}

#define SETBUF_PREAMBLE "csdr"
#define SETBUF_DEFAULT_BUFSIZE 1024
#define STRINGIFY_VALUE(x) STRINGIFY_NAME(x)
#define STRINGIFY_NAME(x) #x

int getbufsize()
{
	if(!env_csdr_dynamic_bufsize_on) return (bigbufs) ? env_csdr_fixed_big_bufsize : env_csdr_fixed_bufsize;
	int recv_first[2];
	fread(recv_first, sizeof(int), 2, stdin);
	if(memcmp(recv_first, SETBUF_PREAMBLE, sizeof(char)*4)!=0)
	{ badsyntax("warning! Did not match preamble on the beginning of the stream. You should put \"csdr setbuf <buffer size>\" at the beginning of the chain! Falling back to default buffer size: " STRINGIFY_VALUE(SETBUF_DEFAULT_BUFSIZE)); return SETBUF_DEFAULT_BUFSIZE; }
	if(recv_first[1]<=0) { badsyntax("warning! Invalid buffer size." ); return 0; }
	return recv_first[1];
}


float* input_buffer;
unsigned char* buffer_u8;
float *output_buffer;
short *buffer_i16;
float *temp_f;
int the_bufsize = 0;
char **argv_global;


#define UNITROUND_UNIT 128

int unitround(int what)
{
	if(what<=0) return UNITROUND_UNIT;
	return ((what-1)&~(UNITROUND_UNIT-1))+UNITROUND_UNIT;
}

int initialize_buffers()
{
	if(!(the_bufsize=getbufsize())) return 0;
	the_bufsize=unitround(the_bufsize);
	if(env_csdr_print_bufsizes) fprintf(stderr,"%s %s: buffer size set to %d\n",argv_global[0], argv_global[1], the_bufsize);
	input_buffer = 	(float*)		malloc(the_bufsize*sizeof(float) * 2); //need the 2× because we might also put complex floats into it
	output_buffer = (float*)		malloc(the_bufsize*sizeof(float) * 2);
	buffer_u8 = 	(unsigned char*)malloc(the_bufsize*sizeof(unsigned char));
	buffer_i16 = 	(short*)		malloc(the_bufsize*sizeof(short));
	temp_f = 		(float*)		malloc(the_bufsize*sizeof(float) * 4);
	return the_bufsize;
}

int sendbufsize(int size)
{
	//The first word is a preamble, "csdr".
	//If the next csdr process detects it, sets the buffer size according to the second word
	if(!env_csdr_dynamic_bufsize_on) return env_csdr_fixed_bufsize;
	if(env_csdr_print_bufsizes) fprintf(stderr,"%s %s: next process proposed input buffer size is %d\n",argv_global[0], argv_global[1], size);
	int send_first[2];
	memcpy((char*)send_first, SETBUF_PREAMBLE, 4*sizeof(char));
	send_first[1] = size;
	fwrite(send_first, sizeof(int), 2, stdout);
	return size;
}

int parse_env()
{
	char* envtmp;
	envtmp=getenv("CSDR_DYNAMIC_BUFSIZE_ON");
	//fprintf(stderr, "envtmp: %s\n",envtmp);
	if(envtmp)
	{
		env_csdr_dynamic_bufsize_on = !!atoi(envtmp);
		env_csdr_fixed_bufsize = 0;
	}
	else
	{
		envtmp=getenv("CSDR_FIXED_BUFSIZE");
		if(envtmp)
		{
			env_csdr_fixed_big_bufsize = env_csdr_fixed_bufsize = atoi(envtmp);
		}
	}
	envtmp=getenv("CSDR_PRINT_BUFSIZES");
	if(envtmp)
	{
		env_csdr_print_bufsizes = atoi(envtmp);
	}
}

int main(int argc, char *argv[])
{
	parse_env();
	argv_global=argv;
	if(argc<=1) return badsyntax(0);
	if(!strcmp(argv[1],"--help")) return badsyntax(0);

	fcntl(STDIN_FILENO, F_SETPIPE_SZ, 65536*32);
	fcntl(STDOUT_FILENO, F_SETPIPE_SZ, 65536*32);
	//fprintf(stderr, "csdr: F_SETPIPE_SZ\n");

	if(!strcmp(argv[1],"setbuf"))
	{
		if(argc<=2) return badsyntax("need required parameter (buffer size)");
		sscanf(argv[2],"%d",&the_bufsize);
		if(the_bufsize<=0) return badsyntax("buffer size <= 0 is invalid");
		sendbufsize(the_bufsize);
		clone_(the_bufsize); //After sending the buffer size out, just copy stdin to stdout
	}

	if(!strcmp(argv[1],"clone"))
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		clone_(the_bufsize);
	}

	if(!strcmp(argv[1],"convert_u8_f"))
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			fread(buffer_u8, sizeof(unsigned char), the_bufsize, stdin);
			convert_u8_f(buffer_u8, output_buffer, the_bufsize);
			FWRITE_R;
			TRY_YIELD;
		}
	}
	if(!strcmp(argv[1],"convert_f_u8")) //not tested
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_R;
			convert_f_u8(input_buffer, buffer_u8, the_bufsize);
			fwrite(buffer_u8, sizeof(unsigned char), the_bufsize, stdout);
			TRY_YIELD;
		}
	}
	if(!strcmp(argv[1],"convert_s8_f"))
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			fread((signed char*)buffer_u8, sizeof(signed char), the_bufsize, stdin);
			convert_s8_f((signed char*)buffer_u8, output_buffer, the_bufsize);
			FWRITE_R;
			TRY_YIELD;
		}
	}
	if(!strcmp(argv[1],"convert_f_s8")) //not tested
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_R;
			convert_f_s8(input_buffer, (signed char*)buffer_u8, the_bufsize);
			fwrite((signed char*)buffer_u8, sizeof(signed char), the_bufsize, stdout);
			TRY_YIELD;
		}
	}
	if((!strcmp(argv[1],"convert_f_i16")) || (!strcmp(argv[1],"convert_f_s16")))
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_R;
			convert_f_i16(input_buffer, buffer_i16, the_bufsize);
			fwrite(buffer_i16, sizeof(short), the_bufsize, stdout);
			TRY_YIELD;
		}
	}
	if((!strcmp(argv[1],"convert_i16_f")) || (!strcmp(argv[1],"convert_s16_f"))) //not tested
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			fread(buffer_i16, sizeof(short), the_bufsize, stdin);
			convert_i16_f(buffer_i16, output_buffer, the_bufsize);
			FWRITE_R;
			TRY_YIELD;
		}
	}
	if(!strcmp(argv[1],"realpart_cf"))
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_C;
			for(int i=0;i<the_bufsize;i++) output_buffer[i]=iof(input_buffer,i);
			FWRITE_R;
			TRY_YIELD;
		}
	}
	if(!strcmp(argv[1],"clipdetect_ff"))
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_R;
			clipdetect_ff(input_buffer, the_bufsize);
			fwrite(input_buffer, sizeof(float), the_bufsize, stdout);
			TRY_YIELD;
		}
	}
	if(!strcmp(argv[1],"gain_ff"))
	{
		if(argc<=2) return badsyntax("need required parameter (gain)");
		float gain;
		sscanf(argv[2],"%g",&gain);
		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_R;
			gain_ff(input_buffer, output_buffer, the_bufsize, gain);
			FWRITE_R;
			TRY_YIELD;
		}
	}
	if(!strcmp(argv[1],"limit_ff"))
	{
		float max_amplitude=1.0;
		if(argc>=3) sscanf(argv[2],"%g",&max_amplitude);
		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_R;
			limit_ff(input_buffer, output_buffer, the_bufsize, max_amplitude);
			FWRITE_R;
			TRY_YIELD;
		}
	}
	if(!strcmp(argv[1],"yes_f"))
	{
		if(argc<=2) return badsyntax("need required parameter (to_repeat)");
		float to_repeat;
		sscanf(argv[2],"%g",&to_repeat);
		int buf_times = 0;
		if(argc>=4) sscanf(argv[3],"%d",&buf_times);
		if(!sendbufsize(initialize_buffers())) return -2;
		for(int i=0;i<the_bufsize;i++) output_buffer[i]=to_repeat;
		for(int i=0;(!buf_times)||i<buf_times;i++)
		{
			fwrite(output_buffer, sizeof(float), the_bufsize, stdout);
			TRY_YIELD;
		}
		return 0;
	}
	if(!strcmp(argv[1],"shift_math_cc"))
	{
		if(argc<=2) return badsyntax("need required parameter (rate)");
		float starting_phase=0;
		float rate;
		sscanf(argv[2],"%g",&rate);
		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			if(!FREAD_C) break;
			FWRITE_C;
			TRY_YIELD;
		}
		return 0;
	}
	//speed tests:
	//csdr yes_f 1 1000000 | time csdr shift_math_cc 0.2 >/dev/null
	//csdr yes_f 1 1000000 | time csdr shift_addition_cc 0.2 >/dev/null
	//csdr yes_f 1 1000000 | time csdr shift_table_cc 0.2 >/dev/null

	if(!strcmp(argv[1],"shift_table_cc"))
	{
		bigbufs=1;
		if(argc<=2) return badsyntax("need required parameter (rate)");
		float starting_phase=0;
		float rate;
		int table_size=65536;
		sscanf(argv[2],"%g",&rate);
		if(argc>3) sscanf(argv[3],"%d",&table_size);
		if(!sendbufsize(initialize_buffers())) return -2;
		shift_table_data_t table_data=shift_table_init(table_size);
		fprintf(stderr,"shift_table_cc: LUT initialized\n");
		for(;;)
		{
			FEOF_CHECK;
			if(!FREAD_C) break;
			starting_phase=shift_table_cc((complexf*)input_buffer, (complexf*)output_buffer, the_bufsize, rate, table_data, starting_phase);
			FWRITE_C;
			TRY_YIELD;
		}
		return 0;
	}

#ifdef LIBCSDR_GPL
	if(!strcmp(argv[1],"decimating_shift_addition_cc"))
	{
		bigbufs=1;
		if(argc<=2) return badsyntax("need required parameter (rate)");
		float starting_phase=0;
		float rate;
		int decimation=1;
		sscanf(argv[2],"%g",&rate);
		if(argc>3) sscanf(argv[3],"%d",&decimation);
		if(!initialize_buffers()) return -2;
		sendbufsize(the_bufsize/decimation);
		shift_addition_data_t d=decimating_shift_addition_init(rate, decimation);
		decimating_shift_addition_status_t s;
		s.decimation_remain=0;
		s.starting_phase=0;
		for(;;)
		{
			FEOF_CHECK;
			if(!FREAD_C) break;
			s=decimating_shift_addition_cc((complexf*)input_buffer, (complexf*)output_buffer, the_bufsize, d, decimation, s);
			fwrite(output_buffer, sizeof(float)*2, s.output_size, stdout);
			TRY_YIELD;
		}
		return 0;
	}

	if(!strcmp(argv[1],"shift_addition_cc"))
	{
		bigbufs=1;

		float starting_phase=0;
		float rate;

		int fd;
		if(fd=init_fifo(argc,argv))
		{
			while(!read_fifo_ctl(fd,"%g\n",&rate)) usleep(10000);
		}
		else
		{
			if(argc<=2) return badsyntax("need required parameter (rate)");
			sscanf(argv[2],"%g",&rate);
		}

		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			shift_addition_data_t data=shift_addition_init(rate);
			fprintf(stderr,"shift_addition_cc: reinitialized to %g\n",rate);
			int remain, current_size;
			float* ibufptr;
			float* obufptr;
			for(;;)
			{
				FEOF_CHECK;
				if(!FREAD_C) break;
				remain=the_bufsize;
				ibufptr=input_buffer;
				obufptr=output_buffer;
				while(remain)
				{
					current_size=(remain>1024)?1024:remain;
					starting_phase=shift_addition_cc((complexf*)ibufptr, (complexf*)obufptr, current_size, data, starting_phase);
					ibufptr+=current_size*2;
					obufptr+=current_size*2;
					remain-=current_size;
				}
				FWRITE_C;
				if(read_fifo_ctl(fd,"%g\n",&rate)) break;
				TRY_YIELD;
			}
		}
		return 0;
	}

	if(!strcmp(argv[1],"shift_addition_cc_test"))
	{
		if(argc<=2) return badsyntax("need required parameter (rate)");
		float rate;
		sscanf(argv[2],"%g",&rate);
		//if(initialize_buffers()) return -2; //most likely we don't need this here
		shift_addition_data_t data=shift_addition_init(rate);
		shift_addition_cc_test(data);
		return 0;
	}
#endif
	if(!strcmp(argv[1],"dcblock_ff"))
	{
		static dcblock_preserve_t dcp; //will be 0 as .bss is set to 0
		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_R;
			dcp=dcblock_ff(input_buffer, output_buffer, the_bufsize, 0, dcp);
			FWRITE_R;
			TRY_YIELD;
		}
	}

	if(!strcmp(argv[1],"fastdcblock_ff"))
	{
		int dcblock_bufsize=SETBUF_DEFAULT_BUFSIZE;
		if(argc>=3) sscanf(argv[2],"%d",&dcblock_bufsize);
		float* dcblock_buffer=(float*)malloc(sizeof(float)*dcblock_bufsize);
		static float last_dc_level=0.0;
		getbufsize(); //it is just dummy
		sendbufsize(dcblock_bufsize);
		for(;;)
		{
			FEOF_CHECK;
			fread(dcblock_buffer, sizeof(float), dcblock_bufsize, stdin);
			last_dc_level=fastdcblock_ff(dcblock_buffer, dcblock_buffer, dcblock_bufsize, last_dc_level);
			fwrite(dcblock_buffer, sizeof(float), dcblock_bufsize, stdout);
			TRY_YIELD;
		}
	}

	if(!strcmp(argv[1],"fmdemod_atan_cf"))
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		float last_phase=0;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_C;
			if(feof(stdin)) return 0;
			last_phase=fmdemod_atan_cf((complexf*)input_buffer, output_buffer, the_bufsize, last_phase);
			FWRITE_R;
			TRY_YIELD;
		}
	}
	if(!strcmp(argv[1],"fmdemod_quadri_cf"))
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		complexf last_sample;
		last_sample.i=0.;
		last_sample.q=0.;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_C;
			last_sample=fmdemod_quadri_cf((complexf*)input_buffer, output_buffer, the_bufsize, temp_f, last_sample);
			FWRITE_R;
			TRY_YIELD;
		}
	}
	if(!strcmp(argv[1],"fmdemod_quadri_novect_cf"))
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		complexf last_sample;
		last_sample.i=0.;
		last_sample.q=0.;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_C;
			last_sample=fmdemod_quadri_novect_cf((complexf*)input_buffer, output_buffer, the_bufsize, last_sample);
			FWRITE_R;
			TRY_YIELD;
		}
	}
	if(!strcmp(argv[1],"deemphasis_wfm_ff"))
	{
		if(argc<=3) return badsyntax("need required parameters (sample rate, tau)");
		if(!sendbufsize(initialize_buffers())) return -2;
		int sample_rate;
		sscanf(argv[2],"%d",&sample_rate);
		float tau;
		sscanf(argv[3],"%g",&tau);
		fprintf(stderr,"deemphasis_wfm_ff: tau = %g, sample_rate = %d\n",tau,sample_rate);
		float last_output=0;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_R;
			last_output=deemphasis_wfm_ff(input_buffer, output_buffer, the_bufsize, tau, sample_rate, last_output);
			FWRITE_R;
			TRY_YIELD;
		}
	}

	if(!strcmp(argv[1],"detect_nan_ff"))
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_R;
			int nan_detect=0;
			for(int i=0; i<the_bufsize;i++)
			{
				if(is_nan(input_buffer[i]))
				{
					nan_detect=1;
					break;
				}
			}
			if(nan_detect) fprintf(stderr, "detect_nan_f: NaN detected!\n");
			fwrite(input_buffer, sizeof(float), the_bufsize, stdout);
			TRY_YIELD;
		}
	}

	if(!strcmp(argv[1],"floatdump_f"))
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_R;
			for(int i=0; i<the_bufsize;i++) fprintf(stderr, "%g ",input_buffer[i]);
			TRY_YIELD;
		}

	}
	if(!strcmp(argv[1],"deemphasis_nfm_ff"))
	{
		if(argc<=2) return badsyntax("need required parameter (sample rate)");
		int sample_rate;
		sscanf(argv[2],"%d",&sample_rate);

		if(!sendbufsize(initialize_buffers())) return -2; //maybe we should take a /2 of bufsize over here

		int processed=0;
		for(;;)
		{
			FEOF_CHECK;
			fread(input_buffer+the_bufsize-processed, sizeof(float), processed, stdin);
			processed=deemphasis_nfm_ff(input_buffer, output_buffer, the_bufsize, sample_rate);
			if(!processed) return badsyntax("deemphasis_nfm_ff: invalid sample rate (this function works only with specific sample rates).");
			memmove(input_buffer,input_buffer+processed,(the_bufsize-processed)*sizeof(float)); //memmove lets the source and destination overlap
			fwrite(output_buffer, sizeof(float), processed, stdout);
			TRY_YIELD;
		}
	}
	if(!strcmp(argv[1],"amdemod_cf"))
	{
		if(!sendbufsize(initialize_buffers())) return -2;

		for(;;)
		{
			FEOF_CHECK;
			FREAD_C;
			amdemod_cf((complexf*)input_buffer, output_buffer, the_bufsize);
			FWRITE_R;
			TRY_YIELD;
		}
	}
	if(!strcmp(argv[1],"amdemod_estimator_cf"))
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_C;
			amdemod_estimator_cf((complexf*)input_buffer, output_buffer, the_bufsize, 0., 0.);
			FWRITE_R;
			TRY_YIELD;
		}
	}
	if(!strcmp(argv[1],"fir_decimate_cc"))
	{
		bigbufs=1;

		if(argc<=2) return badsyntax("need required parameter (decimation factor)");

		int factor;
		sscanf(argv[2],"%d",&factor);

		float transition_bw = 0.05;
		if(argc>=4) sscanf(argv[3],"%g",&transition_bw);

		window_t window = WINDOW_DEFAULT;
		if(argc>=5)
		{
			window=firdes_get_window_from_string(argv[4]);
		}
		else fprintf(stderr,"fir_decimate_cc: window = %s\n",firdes_get_string_from_window(window));

		if(!initialize_buffers()) return -2;
		sendbufsize(the_bufsize/factor);

		int taps_length=firdes_filter_len(transition_bw);
		fprintf(stderr,"fir_decimate_cc: taps_length = %d\n",taps_length);

		int padded_taps_length = taps_length;
		float *taps;
#define NEON_ALIGNMENT (4*4*2)
#ifdef NEON_OPTS
		fprintf(stderr,"taps_length = %d\n", taps_length);
		padded_taps_length = taps_length+(NEON_ALIGNMENT/4)-1 - ((taps_length+(NEON_ALIGNMENT/4)-1)%(NEON_ALIGNMENT/4));
		fprintf(stderr,"padded_taps_length = %d\n", padded_taps_length);

		taps = (float*) (float*)malloc(padded_taps_length+NEON_ALIGNMENT);
		fprintf(stderr,"taps = %x\n", taps);
		taps =  (float*)((((unsigned)taps)+NEON_ALIGNMENT-1) & ~(NEON_ALIGNMENT-1));
		fprintf(stderr,"taps = %x\n", taps);
		for(int i=0;i<padded_taps_length-taps_length;i++) taps[taps_length+i]=0;
#else
		taps=(float*)malloc(taps_length*sizeof(float));
#endif

		firdes_lowpass_f(taps,taps_length,0.5/(float)factor,window);

		int input_skip=0;
		int output_size=0;
		FREAD_C;
		for(;;)
		{
			FEOF_CHECK;
			output_size=fir_decimate_cc((complexf*)input_buffer, (complexf*)output_buffer, the_bufsize, factor, taps, padded_taps_length);
			//fprintf(stderr, "os %d\n",output_size);
			fwrite(output_buffer, sizeof(complexf), output_size, stdout);
			fflush(stdout);
			TRY_YIELD;
			input_skip=factor*output_size;
			memmove((complexf*)input_buffer,((complexf*)input_buffer)+input_skip,(the_bufsize-input_skip)*sizeof(complexf)); //memmove lets the source and destination overlap
			fread(((complexf*)input_buffer)+(the_bufsize-input_skip), sizeof(complexf), input_skip, stdin);
			//fprintf(stderr,"iskip=%d output_size=%d start=%x target=%x skipcount=%x \n",input_skip,output_size,input_buffer, ((complexf*)input_buffer)+(BIG_BUFSIZE-input_skip),(BIG_BUFSIZE-input_skip));
		}
	}
	/*if(!strcmp(argv[1],"ejw_test"))
	{
		printf("ejqd=[");
		complexf ejw;
		float phase=0;
		for(int i=0;i<63;i++)
		{
			e_powj(&ejw,phase);
			phase+=PI*0.3;
			printf("%g+(%g)*i ",iof(&ejw,0),qof(&ejw,0));
		}
		printf("];");
		return 0;
	}*/
	if(!strcmp(argv[1],"firdes_lowpass_f"))
	{
		//Process the params
		if(argc<=3) return badsyntax("need required parameters (cutoff_rate, length)");

		float cutoff_rate;
		sscanf(argv[2],"%g",&cutoff_rate);

		int length;
		sscanf(argv[3],"%d",&length);
		if(length%2==0) return badsyntax("number of symmetric FIR filter taps should be odd");

		window_t window = WINDOW_DEFAULT;
		if(argc>=5)
		{
			window=firdes_get_window_from_string(argv[4]);
		}
		else fprintf(stderr,"firdes_lowpass_f: window = %s\n",firdes_get_string_from_window(window));

		int octave=(argc>=6 && !strcmp("--octave",argv[5]));

		float* taps=(float*)malloc(sizeof(float)*length);

		//Make the filter
		firdes_lowpass_f(taps,length,cutoff_rate,window);

		//Do the output
		if(octave) printf("taps=[");
		for(int i=0;i<length;i++) printf("%g ",taps[i]);
		if(octave) printf("];plot(taps);figure(2);freqz(taps);\n");


		//Wait forever, so that octave won't close just after popping up the window.
		//You can close it with ^C.
		if(octave) { fflush(stdout); getchar(); }
		return 0;
	}
	if(!strcmp(argv[1],"firdes_bandpass_c"))
	{
		//Process the params
		if(argc<=4) return badsyntax("need required parameters (low_cut, high_cut, length)");

		float low_cut;
		sscanf(argv[2],"%g",&low_cut);
		float high_cut;
		sscanf(argv[3],"%g",&high_cut);

		int length;
		sscanf(argv[4],"%d",&length);
		if(length%2==0) return badsyntax("number of symmetric FIR filter taps should be odd");

		window_t window = WINDOW_DEFAULT;
		if(argc>=6)
		{
			window=firdes_get_window_from_string(argv[5]);
		}
		else fprintf(stderr,"firdes_bandpass_c: window = %s\n",firdes_get_string_from_window(window));

		int octave=(argc>=7 && !strcmp("--octave",argv[6]));

		complexf* taps=(complexf*)malloc(sizeof(complexf)*length);

		//Make the filter
		firdes_bandpass_c(taps, length, low_cut, high_cut, window);

		//Do the output
		if(octave) printf("taps=[");
		for(int i=0;i<length;i++) printf("(%g)+(%g)*i ",iof(taps,i),qof(taps,i));
		int fft_length=1024;
		while(fft_length<length) fft_length*=2;
		//if(octave) printf("];\n");
		if(octave) printf(
			"];figure(\"Position\",[0 0 1000 1000]);fser=fft([taps,zeros(1,%d)]);ampl=abs(fser).^2;halfindex=floor(1+size(ampl)(2)/2);\n"
			"amplrev=[ampl(halfindex:end),ampl(1:halfindex)];\n" //we have to swap the output of FFT
			"subplot(2,1,1);plot(amplrev);\n"
			"subplot(2,1,2);plot(arg(fser));\n"
			"#figure(2);freqz(taps);\n"
			"#figur(3);plot3(taps);\n",fft_length-length);

		//Wait forever, so that octave won't close just after popping up the window.
		//You can close it with ^C.
		if(octave) { fflush(stdout); getchar(); }
		return 0;
	}

#ifdef LIBCSDR_GPL
	if(!strcmp(argv[1],"agc_ff"))
	{
		//Process the params
		//Explanation of what these actually do is in the DSP source.
		//These good default values are for SSB sampled at 48000 kHz.
		short hang_time=200;
		if(argc>=3) sscanf(argv[2],"%hd",&hang_time);

		float reference=0.2;
		if(argc>=4) sscanf(argv[3],"%g",&reference);

		float attack_rate=0.01;
		if(argc>=5) sscanf(argv[4],"%g",&attack_rate);

		float decay_rate=0.0001;
		if(argc>=6) sscanf(argv[5],"%g",&decay_rate);

		float max_gain=65536;
		if(argc>=7) sscanf(argv[6],"%g",&max_gain);

		short attack_wait=0;
		if(argc>=8) sscanf(argv[7],"%hd",&attack_wait);

		float filter_alpha=0.999;//0.001;
		if(argc>=9) sscanf(argv[8],"%g",&filter_alpha);

		if(!sendbufsize(initialize_buffers())) return -2;

		float last_gain=1.0;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_R;
			last_gain=agc_ff(input_buffer, output_buffer, the_bufsize, reference, attack_rate, decay_rate, max_gain, hang_time, attack_wait, filter_alpha, last_gain);
			FWRITE_R;
			TRY_YIELD;
		}
	}
#endif

	if(!strcmp(argv[1],"fastagc_ff"))
	{

		static fastagc_ff_t input; //is in .bss and gets cleared to zero before main()

		input.input_size=1024;
		if(argc>=3) sscanf(argv[2],"%d",&input.input_size);

		getbufsize(); //dummy
		sendbufsize(input.input_size);

		input.reference=1.0;
		if(argc>=4) sscanf(argv[3],"%g",&input.reference);

		//input.max_peak_ratio=12.0;
		//if(argc>=5) sscanf(argv[3],"%g",&input.max_peak_ratio);

		input.buffer_1=(float*)calloc(input.input_size,sizeof(float));
		input.buffer_2=(float*)calloc(input.input_size,sizeof(float));
		input.buffer_input=(float*)malloc(sizeof(float)*input.input_size);
		float* agc_output_buffer=(float*)malloc(sizeof(float)*input.input_size);
		for(;;)
		{
			FEOF_CHECK;
			fread(input.buffer_input, sizeof(float), input.input_size, stdin);
			fastagc_ff(&input, agc_output_buffer);
			fwrite(agc_output_buffer, sizeof(float), input.input_size, stdout);
			TRY_YIELD;
		}
	}

	int suboptimal;
	if( (suboptimal=!strcmp(argv[1],"suboptimal_rational_resampler_ff"))||(!strcmp(argv[1],"rational_resampler_ff")) )
	{

		//last@2014-11-06: ./docompile; ./csdr yes_f 1.0 | ./csdr suboptimal_rational_resampler_ff 5 2

		//Process the params
		if(argc<=3) return badsyntax("need required parameters (interpolation, decimation)");
		int interpolation;
		sscanf(argv[2],"%d",&interpolation);
		int decimation;
		sscanf(argv[3],"%d",&decimation);

		float transition_bw=0.05;
		if(argc>=5) sscanf(argv[4],"%g",&transition_bw);

		window_t window = WINDOW_DEFAULT;
		if(argc>=6)
		{
			window=firdes_get_window_from_string(argv[5]);
		}
		else fprintf(stderr,"rational_resampler_ff: window = %s\n",firdes_get_string_from_window(window));

		if(suboptimal) fprintf(stderr,"note: suboptimal rational resampler chosen.\n");

		if(!initialize_buffers()) return -2;

		if(decimation==1&&interpolation==1) { sendbufsize(the_bufsize); clone_(the_bufsize); } //copy input to output in this special case (and stick in this function).

		//Alloc output buffer
		int resampler_output_buffer_size=(the_bufsize*interpolation)/decimation;
		sendbufsize(resampler_output_buffer_size);
		float* resampler_output_buffer=(float*)malloc(sizeof(float)*resampler_output_buffer_size);
		float* suboptimal_resampler_temp_buffer = (suboptimal)?(float*)malloc(sizeof(float)*the_bufsize*interpolation):NULL;

		//Generate filter taps
		int taps_length = firdes_filter_len(transition_bw);
		float* taps = (float*)malloc(sizeof(float)*taps_length);
		rational_resampler_get_lowpass_f(taps, taps_length, interpolation, decimation, window);

		static rational_resampler_ff_t d; //in .bss => initialized to zero

		for(;;)
		{
			FEOF_CHECK;
			if(d.input_processed==0) d.input_processed=the_bufsize;
			else memcpy(input_buffer, input_buffer+d.input_processed, sizeof(float)*(the_bufsize-d.input_processed));
			fread(input_buffer+(the_bufsize-d.input_processed), sizeof(float), d.input_processed, stdin);
			//if(suboptimal) d=suboptimal_rational_resampler_ff(input_buffer, resampler_output_buffer, the_bufsize, interpolation, decimation, taps, taps_length, suboptimal_resampler_temp_buffer); else
			d=rational_resampler_ff(input_buffer, resampler_output_buffer, the_bufsize, interpolation, decimation, taps, taps_length, d.last_taps_delay);
			//fprintf(stderr,"resampled %d %d, %d\n",d.output_size, d.input_processed, d.input_processed);
			fwrite(resampler_output_buffer, sizeof(float), d.output_size, stdout);
			TRY_YIELD;
		}
	}


	if(!strcmp(argv[1],"fractional_decimator_ff"))
	{
		//Process the params
		if(argc<=2) return badsyntax("need required parameters (rate)");
		float rate;
		sscanf(argv[2],"%g",&rate);

		float transition_bw=0.03;
		if(argc>=4) sscanf(argv[3],"%g",&transition_bw);

		window_t window = WINDOW_DEFAULT;
		if(argc>=5)
		{
			window = firdes_get_window_from_string(argv[4]);
		}
		else fprintf(stderr,"fractional_decimator_ff: window = %s\n",firdes_get_string_from_window(window));

		if(!initialize_buffers()) return -2;
		sendbufsize(the_bufsize / rate);

		if(rate==1) clone_(the_bufsize); //copy input to output in this special case (and stick in this function).

		//Generate filter taps
		int taps_length = firdes_filter_len(transition_bw);
		fprintf(stderr,"fractional_decimator_ff: taps_length = %d\n",taps_length);
		float* taps = (float*)malloc(sizeof(float)*taps_length);
		firdes_lowpass_f(taps, taps_length, 0.59*0.5/(rate-transition_bw), window); //0.6 const to compensate rolloff
		//for(int=0;i<taps_length; i++) fprintf(stderr,"%g ",taps[i]);

		static fractional_decimator_ff_t d; //in .bss => initialized to zero
		for(;;)
		{
			FEOF_CHECK;
			if(d.input_processed==0) d.input_processed=the_bufsize;
			else memcpy(input_buffer, input_buffer+d.input_processed, sizeof(float)*(the_bufsize-d.input_processed));
			fread(input_buffer+(the_bufsize-d.input_processed), sizeof(float), d.input_processed, stdin);
			d = fractional_decimator_ff(input_buffer, output_buffer, the_bufsize, rate, taps, taps_length, d);
			fwrite(output_buffer, sizeof(float), d.output_size, stdout);
			TRY_YIELD;
		}
	}

	if(!strcmp(argv[1],"fft_cc"))
	{
		if(argc<=3) return badsyntax("need required parameters (fft_size, out_of_every_n_samples)");
		int fft_size;
		sscanf(argv[2],"%d",&fft_size);
		if(log2n(fft_size)==-1) return badsyntax("fft_size should be power of 2");
		int every_n_samples;
		sscanf(argv[3],"%d",&every_n_samples);
		int benchmark=0;
		int octave=0;
		window_t window = WINDOW_DEFAULT;
		if(argc>=5)
		{
			window=firdes_get_window_from_string(argv[4]);
		}
		if(argc>=6)
		{
			benchmark|=!strcmp("--benchmark",argv[5]);
			octave|=!strcmp("--octave",argv[5]);
		}
		if(argc>=7)
		{
			benchmark|=!strcmp("--benchmark",argv[6]);
			octave|=!strcmp("--octave",argv[6]);
		}

		if(!initialize_buffers()) return -2;
		sendbufsize(fft_size);

		//make FFT plan
		complexf* input=(complexf*)fft_malloc(sizeof(complexf)*fft_size);
		complexf* windowed=(complexf*)fft_malloc(sizeof(complexf)*fft_size);
		complexf* output=(complexf*)fft_malloc(sizeof(complexf)*fft_size);
		if(benchmark) fprintf(stderr,"fft_cc: benchmarking...");
		FFT_PLAN_T* plan=make_fft_c2c(fft_size, windowed, output, 1, benchmark);
		if(benchmark) fprintf(stderr," done\n");
		if(octave) printf("setenv(\"GNUTERM\",\"X11 noraise\");y=zeros(1,%d);semilogy(y,\"ydatasource\",\"y\");\n",fft_size);
		for(;;)
		{
			FEOF_CHECK;
			if(every_n_samples>fft_size)
			{
				fread(input, sizeof(complexf), fft_size, stdin);
				//skipping samples before next FFT (but fseek doesn't work for pipes)
				for(int seek_remain=every_n_samples-fft_size;seek_remain>0;seek_remain-=the_bufsize)
				{
					fread(temp_f, sizeof(complexf), MIN_M(the_bufsize,seek_remain), stdin);
				}
			}
			else
			{
				//overlapped FFT
				for(int i=0;i<fft_size-every_n_samples;i++) input[i]=input[i+every_n_samples];
				fread(input+fft_size-every_n_samples, sizeof(complexf), every_n_samples, stdin);
			}
			apply_window_c(input,windowed,fft_size,window);
			fft_execute(plan);
			if(octave)
			{
				printf("fftdata=[");
				//we have to swap the two parts of the array to get a valid spectrum
				for(int i=fft_size/2;i<fft_size;i++) printf("(%g)+(%g)*i ",iof(output,i),qof(output,i));
				for(int i=0;i<fft_size/2;i++) printf("(%g)+(%g)*i ",iof(output,i),qof(output,i));
				printf(
					"];\n"
					"y=abs(fftdata);\n"
					"refreshdata;\n"
				);
			}
			else fwrite(output, sizeof(complexf), fft_size, stdout);
			TRY_YIELD;
		}
	}
	#define LOGPOWERCF_BUFSIZE 64
	if(!strcmp(argv[1],"logpower_cf"))
	{
		float add_db=0;
		if(argc>=3) sscanf(argv[2],"%g",&add_db);

		if(!sendbufsize(initialize_buffers())) return -2;

		for(;;)
		{
			FEOF_CHECK;
			fread(input_buffer, sizeof(complexf), the_bufsize, stdin);
			logpower_cf((complexf*)input_buffer,output_buffer, the_bufsize, add_db);
			fwrite(output_buffer, sizeof(float), the_bufsize, stdout);
			TRY_YIELD;
		}
	}

	if(!strcmp(argv[1],"fft_exchange_sides_ff"))
	{
		if(argc<=2) return badsyntax("need required parameters (fft_size)");
		int fft_size;
		sscanf(argv[2],"%d",&fft_size);
		if(!getbufsize()) return -2; //dummy
		sendbufsize(fft_size);
		float* input_buffer_s1 = (float*)malloc(sizeof(float)*fft_size/2);
		float* input_buffer_s2 = (float*)malloc(sizeof(float)*fft_size/2);
		for(;;)
		{
			FEOF_CHECK;
			fread(input_buffer_s1, sizeof(float), fft_size/2, stdin);
			fread(input_buffer_s2, sizeof(float), fft_size/2, stdin);
			fwrite(input_buffer_s2, sizeof(float), fft_size/2, stdout);
			fwrite(input_buffer_s1, sizeof(float), fft_size/2, stdout);
			TRY_YIELD;
		}
	}


#ifdef USE_IMA_ADPCM

#define COMPRESS_FFT_PAD_N 10
//We will pad the FFT at the beginning, with the first value of the input data, COMPRESS_FFT_PAD_N times.
//No, this is not advanced DSP, just the ADPCM codec produces some gabarge samples at the beginning,
//so we just add data to become garbage and get skipped.
//COMPRESS_FFT_PAD_N should be even.

	if(!strcmp(argv[1],"compress_fft_adpcm_f_u8"))
	{
		if(argc<=2) return badsyntax("need required parameters (fft_size)");
		int fft_size;
		sscanf(argv[2],"%d",&fft_size);
		int real_data_size=fft_size+COMPRESS_FFT_PAD_N;
		if(!getbufsize()) return -2; //dummy
		sendbufsize(real_data_size);
		float* input_buffer_cwa = (float*)malloc(sizeof(float)*real_data_size);
		short* temp_buffer_cwa = (short*)malloc(sizeof(short)*real_data_size);
		unsigned char* output_buffer_cwa = (unsigned char*)malloc(sizeof(unsigned char)*(real_data_size/2));
		ima_adpcm_state_t d;
		d.index=d.previousValue=0;
		for(;;)
		{
			FEOF_CHECK;
			fread(input_buffer_cwa+COMPRESS_FFT_PAD_N, sizeof(float), fft_size, stdin);
			for(int i=0;i<COMPRESS_FFT_PAD_N;i++) input_buffer_cwa[i]=input_buffer_cwa[COMPRESS_FFT_PAD_N]; //do padding
			for(int i=0;i<real_data_size;i++) temp_buffer_cwa[i]=input_buffer_cwa[i]*100; //convert float dB values to short
			encode_ima_adpcm_i16_u8(temp_buffer_cwa, output_buffer_cwa, real_data_size, d); //we always return to original d at any new buffer
			fwrite(output_buffer_cwa, sizeof(unsigned char), real_data_size/2, stdout);
			TRY_YIELD;
		}
	}
#endif

#define TIME_TAKEN(start,end) ((end.tv_sec-start.tv_sec)+(end.tv_nsec-start.tv_nsec)/1e9)

	if(!strcmp(argv[1],"fft_benchmark"))
	{
		if(argc<=3) return badsyntax("need required parameters (fft_size, fft_cycles)");
		int fft_size;
		sscanf(argv[2],"%d",&fft_size);
		int fft_cycles;
		sscanf(argv[3],"%d",&fft_cycles);

		int benchmark=(argc>=5)&&!strcmp(argv[4],"--benchmark");
		fprintf(stderr,"fft_benchmark: FFT library used: %s\n",FFT_LIBRARY_USED);

		complexf* input=(complexf*)fft_malloc(sizeof(complexf)*fft_size);
		complexf* output=(complexf*)fft_malloc(sizeof(complexf)*fft_size);

		//fill input with random data
		srand(time(NULL));
		for(int i=0;i<fft_size;i++)
		{
			iof(input,i)=rand()/(float)INT_MAX;
			qof(input,i)=rand()/(float)INT_MAX;
		}

		//initialize FFT library, and measure time
		fprintf(stderr,"fft_benchmark: initializing... ");
		struct timespec start_time, end_time;
		clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
		FFT_PLAN_T* plan=make_fft_c2c(fft_size,input,output,1,benchmark);
		clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);
		fprintf(stderr,"done in %g seconds.\n",TIME_TAKEN(start_time,end_time));

		//do the actual measurement about the FFT
		clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
		for(int i=0;i<fft_cycles;i++) fft_execute(plan);
		clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);
		float time_taken_fft = TIME_TAKEN(start_time,end_time);
		fprintf(stderr,"fft_benchmark: %d transforms of %d processed in %g seconds, %g seconds each.\n",fft_cycles,fft_size,time_taken_fft,time_taken_fft/fft_cycles);
		return 0;
	}

	if(!strcmp(argv[1],"bandpass_fir_fft_cc")) //this command does not exist as a separate function
	{
		float low_cut;
		float high_cut;
		float transition_bw;
		window_t window = WINDOW_DEFAULT;
		char window_string[256]; //TODO: nice buffer overflow opportunity

		int fd;
		if(fd=init_fifo(argc,argv))
		{
			while(!read_fifo_ctl(fd,"%g %g\n",&low_cut,&high_cut)) usleep(10000);
			if(argc<=4) return badsyntax("need more required parameters (transition_bw)");
		}
		else
		{
			if(argc<=4) return badsyntax("need required parameters (low_cut, high_cut, transition_bw)");
			sscanf(argv[2],"%g",&low_cut);
			sscanf(argv[3],"%g",&high_cut);
		}
		sscanf(argv[4],"%g",&transition_bw);
		if(argc>=6)	window=firdes_get_window_from_string(argv[5]);
		else fprintf(stderr,"bandpass_fir_fft_cc: window = %s\n",firdes_get_string_from_window(window));

		//calculate the FFT size and the other length parameters
		int taps_length=firdes_filter_len(transition_bw); //the number of non-zero taps
		int fft_size=next_pow2(taps_length); //we will have to pad the taps with zeros until the next power of 2 for FFT
		//the number of padding zeros is the number of output samples we will be able to take away after every processing step, and it looks sane to check if it is large enough.
		if (fft_size-taps_length<200) fft_size<<=1;
		int input_size = fft_size - taps_length + 1;
		int overlap_length = taps_length - 1;
		fprintf(stderr,"bandpass_fir_fft_cc: (fft_size = %d) = (taps_length = %d) + (input_size = %d) - 1\n(overlap_length = %d) = taps_length - 1\n", fft_size, taps_length, input_size, overlap_length);
		if (fft_size<=2) return badsyntax("FFT size error.");

		if(!sendbufsize(getbufsize())) return -2;

		//prepare making the filter and doing FFT on it
		complexf* taps=(complexf*)calloc(sizeof(complexf),fft_size); //initialize to zero
		complexf* taps_fft=(complexf*)malloc(sizeof(complexf)*fft_size);
		FFT_PLAN_T* plan_taps = make_fft_c2c(fft_size, taps, taps_fft, 1, 0); //forward, don't benchmark (we need this only once)

		//make FFT plans for continously processing the input
		complexf* input = fft_malloc(fft_size*sizeof(complexf));
		complexf* input_fourier = fft_malloc(fft_size*sizeof(complexf));
		FFT_PLAN_T* plan_forward = make_fft_c2c(fft_size, input, input_fourier, 1, 1); //forward, do benchmark

		complexf* output_fourier = fft_malloc(fft_size*sizeof(complexf));
		complexf* output_1 = fft_malloc(fft_size*sizeof(complexf));
		complexf* output_2 = fft_malloc(fft_size*sizeof(complexf));
		//we create 2x output buffers so that one will preserve the previous overlap:
		FFT_PLAN_T* plan_inverse_1 = make_fft_c2c(fft_size, output_fourier, output_1, 0, 1); //inverse, do benchmark
		FFT_PLAN_T* plan_inverse_2 = make_fft_c2c(fft_size, output_fourier, output_2, 0, 1);
		//we initialize this buffer to 0 as it will be taken as the overlap source for the first time:
		for(int i=0;i<fft_size;i++) iof(plan_inverse_2->output,i)=qof(plan_inverse_2->output,i)=0;

		for(int i=input_size;i<fft_size;i++) iof(input,i)=qof(input,i)=0; //we pre-pad the input buffer with zeros

		for(;;)
		{
			//make the filter
			fprintf(stderr,"bandpass_fir_fft_cc: filter initialized, low_cut = %g, high_cut = %g\n",low_cut,high_cut);
			firdes_bandpass_c(taps, taps_length, low_cut, high_cut, window);
			fft_execute(plan_taps);

			for(int odd=0;;odd=!odd) //the processing loop
			{
				FEOF_CHECK;
				fread(input, sizeof(complexf), input_size, stdin);
				FFT_PLAN_T* plan_inverse = (odd)?plan_inverse_2:plan_inverse_1;
				FFT_PLAN_T* plan_contains_last_overlap = (odd)?plan_inverse_1:plan_inverse_2; //the other
				complexf* last_overlap = (complexf*)plan_contains_last_overlap->output + input_size; //+ fft_size - overlap_length;
				apply_fir_fft_cc (plan_forward, plan_inverse, taps_fft, last_overlap, overlap_length);
				int returned=fwrite(plan_inverse->output, sizeof(complexf), input_size, stdout);
				if(read_fifo_ctl(fd,"%g %g\n",&low_cut,&high_cut)) break;
				TRY_YIELD;
			}
		}

	}

#ifdef USE_IMA_ADPCM
#define IMA_ADPCM_BUFSIZE BUFSIZE

	if( (!strcmp(argv[1],"encode_ima_adpcm_i16_u8"))||(!strcmp(argv[1],"encode_ima_adpcm_s16_u8")) )
	{
		if(!sendbufsize(initialize_buffers()/2)) return -2;
		ima_adpcm_state_t d;
		d.index=d.previousValue=0;
		for(;;)
		{
			FEOF_CHECK;
			fread(buffer_i16, sizeof(short), the_bufsize, stdin);
			d=encode_ima_adpcm_i16_u8(buffer_i16, buffer_u8, the_bufsize, d);
			fwrite(buffer_u8, sizeof(unsigned char), the_bufsize/2, stdout);
			TRY_YIELD;
		}
	}

	if( (!strcmp(argv[1],"decode_ima_adpcm_u8_i16"))||(!strcmp(argv[1],"decode_ima_adpcm_u8_s16")) )
	{
		ima_adpcm_state_t d;
		d.index=d.previousValue=0;
		if(!sendbufsize(initialize_buffers()*2)) return -2;
		for(;;)
		{
			FEOF_CHECK;
			fread(buffer_u8, sizeof(unsigned char), the_bufsize, stdin);
			d=decode_ima_adpcm_u8_i16(buffer_u8, buffer_i16, the_bufsize, d);
			fwrite(buffer_i16, sizeof(short), the_bufsize*2, stdout);
			TRY_YIELD;
		}
	}
#endif

	if(!strcmp(argv[1],"flowcontrol"))
	{
		if(argc<=3) return badsyntax("need required parameters (data_rate, reads_per_seconds)");
		int data_rate;
		sscanf(argv[2],"%d",&data_rate);
		int reads_per_second;
		sscanf(argv[3],"%d",&reads_per_second);
		int flowcontrol_bufsize=ceil(1.*(double)data_rate/reads_per_second);
		if(!getbufsize()) return -2;
		sendbufsize(flowcontrol_bufsize);
		unsigned char* flowcontrol_buffer = (unsigned char*)malloc(sizeof(unsigned char)*flowcontrol_bufsize);
		int flowcontrol_sleep=floor(1000000./reads_per_second);
		fprintf(stderr, "flowcontrol: flowcontrol_bufsize = %d, flowcontrol_sleep = %d\n", flowcontrol_bufsize, flowcontrol_sleep);
		for(;;)
		{
			FEOF_CHECK;
			fread(flowcontrol_buffer, sizeof(unsigned char), flowcontrol_bufsize, stdin);
			fwrite(flowcontrol_buffer, sizeof(unsigned char), flowcontrol_bufsize, stdout);
			usleep(flowcontrol_sleep);
			TRY_YIELD;
		}
	}

#if 0
	if(!strcmp(argv[1],"flowcontrol"))
	{
		if(argc<=3) return badsyntax("need required parameters (data_rate, reads_per_seconds)");

		int data_rate;
		sscanf(argv[2],"%d",&data_rate);

		int reads_per_second=0;
		if(strcmp(argv[3],"auto")) sscanf(argv[3],"%d",&reads_per_second);

		float prebuffer=2;
		if(argc>4) sscanf(argv[4],"%g",&prebuffer);

		int thrust=10;
		if(argc>5) sscanf(argv[5],"%d",&thrust);

		int flowcontrol_readsize, flowcontrol_bufsize, got_bufsize;

		if(!(got_bufsize=getbufsize())) return -2;

		if(reads_per_second)
		{
			flowcontrol_readsize=ceil(1.*(double)data_rate/reads_per_second);
		}
		else
		{
			flowcontrol_readsize=got_bufsize;
			reads_per_second=data_rate/flowcontrol_readsize;
		}
		flowcontrol_bufsize=flowcontrol_readsize*floor(reads_per_second*prebuffer);

		int flowcontrol_bufindex=0;
		unsigned char* flowcontrol_buffer = (unsigned char*)malloc(sizeof(unsigned char)*flowcontrol_bufsize);
		int flowcontrol_sleep=floor(1000000./reads_per_second);

		fcntl(STDIN_FILENO, F_SETFL, fcntl(STDIN_FILENO, F_GETFL, 0) | O_NONBLOCK);

		sendbufsize(flowcontrol_readsize);
		fflush(stdout);

		int flowcontrol_is_buffering = 1;
		int read_return;

		struct timespec start_time, end_time;

		unsigned long long int all_bytes_written=0;
		int test=0;

		fprintf(stderr, "flowcontrol: flowcontrol_readsize = %d, flowcontrol_bufsize = %d, flowcontrol_sleep = %d\n", flowcontrol_readsize, flowcontrol_bufsize, flowcontrol_sleep);
		for (; ;) //my friend has told me that this is like two smileys ;)
		{
			FEOF_CHECK;
			fprintf(stderr, "r");
			read_return=read(STDIN_FILENO, flowcontrol_buffer+flowcontrol_bufindex, sizeof(unsigned char) * (flowcontrol_bufsize-flowcontrol_bufindex) );
			fprintf(stderr, "t");
			if(read_return>0) flowcontrol_bufindex+=read_return;


			if(flowcontrol_is_buffering)
			{
				fprintf(stderr, "flowcontrol: buffering, flowcontrol_bufindex = %d\n", flowcontrol_bufindex);
				if(flowcontrol_bufindex==flowcontrol_bufsize) { flowcontrol_is_buffering = 0; clock_gettime(CLOCK_MONOTONIC_RAW, &start_time); }
				else if(read_return<=0) continue;
			}
			else {
				clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);
				int thrust_added=0;
				while( (all_bytes_written+thrust*flowcontrol_readsize) / TIME_TAKEN(start_time,end_time) < data_rate )
				{
					thrust_added |= thrust++;
				}
				//if(!(test++%10)) fprintf(stderr, "abw=%g\n", all_bytes_written / TIME_TAKEN(start_time,end_time));
				/*if(!thrust_added && TIME_TAKEN(start_time,end_time)>50)
				{
					clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
					all_bytes_written=0;
				}*/
				while(all_bytes_written>data_rate && TIME_TAKEN(start_time,end_time)>1)
				{
					all_bytes_written-=data_rate;
					start_time.tv_sec++;
				}
				do
				{
					//if(thrust) fprintf(stderr, "flowcontrol: %d .. thrust\n", thrust);
					write(STDOUT_FILENO, flowcontrol_buffer, flowcontrol_readsize);
					fflush(stdout);
					//fsync(STDOUT_FILENO);
					memmove(flowcontrol_buffer, flowcontrol_buffer+flowcontrol_readsize, flowcontrol_bufindex-flowcontrol_readsize);
					flowcontrol_bufindex -= flowcontrol_readsize;
					all_bytes_written += flowcontrol_readsize;
				} while(thrust && thrust-- && flowcontrol_bufindex>=flowcontrol_readsize);
			}

			usleep(flowcontrol_sleep);
			TRY_YIELD;
		}
	}
#endif

	if(!strcmp(argv[1],"through"))
	{
		struct timespec start_time, end_time;
		if(!sendbufsize(initialize_buffers())) return -2;

		int time_now_sec=0;
		int buffer_count=0;

		unsigned char* through_buffer;
		through_buffer = (unsigned char*)malloc(the_bufsize*sizeof(float));


		for(;;)
		{
			FEOF_CHECK;
			fread(through_buffer, sizeof(float), the_bufsize, stdin);

			if(!time_now_sec)
			{
				time_now_sec=1;
				clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
			}
			else
			{
				clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);
				float timetaken;
				if(time_now_sec<(timetaken=TIME_TAKEN(start_time,end_time)))
				{
					fprintf( stderr, "through: %lu bytes/s %d\n", (unsigned long)floor((float)buffer_count*the_bufsize*sizeof(float)/timetaken), buffer_count );
					time_now_sec=ceil(timetaken);
				}
			}
			fwrite(through_buffer, sizeof(float), the_bufsize, stdout);
			buffer_count++;
			TRY_YIELD;
		}
	}

	if(!strcmp(argv[1],"dsb_fc"))
	{
		float q_value = 0;
		if(argc>=3) sscanf(argv[2],"%g",&q_value);

		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_R;
			for(int i=0;i<the_bufsize;i++)
			{
				iof(output_buffer,i)=input_buffer[i];
				qof(output_buffer,i)=q_value;
			}
			FWRITE_C;
			TRY_YIELD;
		}
	}

	if(!strcmp(argv[1],"convert_f_samplerf"))
	{
		if(argc<=2) return badsyntax("need required parameter (wait_for_this_sample)");

		unsigned wait_for_this_sample;
		sscanf(argv[2],"%u",&wait_for_this_sample);

		if(!sendbufsize(initialize_buffers())) return -2;
		unsigned char* samplerf_buf = (unsigned char*) malloc(16*the_bufsize);
		for(;;)
		{
			FEOF_CHECK;
			FREAD_R;
			for(int i=0;i<the_bufsize;i++)
			{
				*((double*)(&samplerf_buf[16*i])) = input_buffer[i];
				*((unsigned*)(&samplerf_buf[16*i+8])) = wait_for_this_sample;
				*((unsigned*)(&samplerf_buf[16*i+12])) = 0;

			}
			fwrite(samplerf_buf, 16, the_bufsize, stdout);
			TRY_YIELD;
		}
	}

	if(!strcmp(argv[1],"add_dcoffset_cc"))
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_C;
			add_dcoffset_cc((complexf*)input_buffer, (complexf*)output_buffer, the_bufsize);
			FWRITE_C;
			TRY_YIELD;
		}
	}

	if(!strcmp(argv[1],"fmmod_fc"))
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		float last_phase = 0;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_R;
			last_phase = fmmod_fc(input_buffer, (complexf*)output_buffer, the_bufsize, last_phase);
			FWRITE_C;
			TRY_YIELD;
		}
	}

	if(!strcmp(argv[1],"fixed_amplitude_cc"))
	{
		if(argc<=2) return badsyntax("need required parameter (new_amplitude)");

		float new_amplitude;
		sscanf(argv[2],"%g",&new_amplitude);

		if(!sendbufsize(initialize_buffers())) return -2;
		for(;;)
		{
			FEOF_CHECK;
			FREAD_C;
			fixed_amplitude_cc((complexf*)input_buffer, (complexf*)output_buffer, the_bufsize, new_amplitude);
			FWRITE_C;
			TRY_YIELD;
		}
	}

	if((!strcmp(argv[1],"mono2stereo_i16"))||(!strcmp(argv[1],"mono2stereo_s16")))
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		float last_phase = 0;
		for(;;)
		{
			FEOF_CHECK;
			fread (input_buffer, sizeof(short), the_bufsize, stdin);
			for(int i=0;i<the_bufsize;i++)
			{
				*(((short*)output_buffer)+2*i)=*(((short*)input_buffer)+i);
				*(((short*)output_buffer)+2*i+1)=*(((short*)input_buffer)+i);
			}
			fwrite (output_buffer, sizeof(short)*2, the_bufsize, stdout);
			TRY_YIELD;
		}
	}

	if(!strcmp(argv[1],"squelch_and_smeter_cc"))
	{
		if(!sendbufsize(initialize_buffers())) return -2;
		float power;
		float squelch_level;
		int decimation;
		int report_every_nth;
		int fd;
		char power_value_buf[101];
		int power_value_buf_size;
		int report_cntr=0;
		complexf* zerobuf = (complexf*)malloc(sizeof(complexf)*the_bufsize);
		for(int i=0;i<the_bufsize*2;i++) *(((float*)zerobuf)+i)=0;
		if(fd=init_fifo(argc,argv)) while(!read_fifo_ctl(fd,"%g\n",&squelch_level)) usleep(10000);
		else return badsyntax("need required parameter (--fifo <fifo>)");
		fprintf(stderr, "squelch_and_power_cc: initial squelch level is %g\n", squelch_level);
		if((argc<=5)||((argc>5)&&(strcmp(argv[4],"--outfifo")))) return badsyntax("need required parameter (--outfifo <fifo>)");
		int fd2 = open(argv[5], O_WRONLY);
		if(fd2==-1) return badsyntax("error while opening --outfifo");
		int flags = fcntl(fd2, F_GETFL, 0);
		fcntl(fd2, F_SETFL, flags | O_NONBLOCK);
		if(argc<=6) return badsyntax("need required parameter (use_every_nth)");
		sscanf(argv[6],"%d",&decimation);
		if(decimation<=0) return badsyntax("use_every_nth <= 0 is invalid");
		if(argc<=7) return badsyntax("need required parameter (report_every_nth)");
		sscanf(argv[7],"%d",&report_every_nth);
		if(report_every_nth<=0) return badsyntax("report_every_nth <= 0 is invalid");
		for(;;)
		{
			FEOF_CHECK;
			FREAD_C; //read input data
			power = get_power_c((complexf*)input_buffer, the_bufsize, decimation);
			if(report_cntr++>report_every_nth)
			{
				report_cntr=0;
				power_value_buf_size=snprintf(power_value_buf,100,"%g\n",power);
				write(fd2,power_value_buf,power_value_buf_size*sizeof(char));
		  }
			if(squelch_level==0||power>=squelch_level)
			{
				//fprintf(stderr,"P");
				fwrite(input_buffer, sizeof(complexf), the_bufsize, stdout);
			}
			else
			{
				//fprintf(stderr,"S");
				fwrite(zerobuf, sizeof(complexf), the_bufsize, stdout);
			}
			if(read_fifo_ctl(fd,"%g\n",&squelch_level)) fprintf(stderr, "squelch_and_power_cc: new squelch level is %g\n", squelch_level);
			TRY_YIELD;
		}
	}

	if(!strcmp(argv[1],"none"))
	{
		return 0;
	}

	return badsyntax("function name given in argument 1 does not exist. Possible causes:\n- You mistyped the commandline.\n- You need to update csdr to a newer version (if available).");

}
