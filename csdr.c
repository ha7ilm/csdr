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
#include <strings.h>
#include <errno.h>
#include "fastddc.h"
#include <assert.h>

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
"    convert_f_s24 [--bigendian]\n"
"    convert_s24_f [--bigendian]\n"
"    realpart_cf\n"
"    clipdetect_ff\n"
"    limit_ff [max_amplitude]\n"
"    gain_ff <gain>\n"
"    clone\n"
"    none\n"
"    yes_f <to_repeat> [buf_times]\n"
"    detect_nan_ff\n"
"    dump_f\n"
"    shift_math_cc <rate>\n"
"    shift_math_cc --fifo <fifo_path>\n"
"    shift_addition_cc <rate>\n"
"    shift_addition_cc --fifo <fifo_path>\n"
"    shift_addition_cc_test\n"
"    shift_table_cc <rate> [table_size]\n"
"    shift_addition_fc <rate>\n"
"    shift_addition_fc --fifo <fifo_path>\n"
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
"    fir_interpolate_cc <interpolation_factor> [transition_bw [window]]\n"
"    firdes_lowpass_f <cutoff_rate> <length> [window [--octave]]\n"
"    firdes_bandpass_c <low_cut> <high_cut> <length> [window [--octave]]\n"
"    agc_ff [hang_time [reference [attack_rate [decay_rate [max_gain [attack_wait [filter_alpha]]]]]]]\n"
"    fastagc_ff [block_size [reference]]\n"
"    rational_resampler_ff <interpolation> <decimation> [transition_bw [window]]\n"
"    old_fractional_decimator_ff <decimation_rate> [transition_bw [window]]\n"
"    fractional_decimator_ff <decimation_rate> [num_poly_points ( [transition_bw [window]] | --prefilter )]\n"
"    fft_cc <fft_size> <out_of_every_n_samples> [window [--octave] [--benchmark]]\n"
"    fft_fc <fft_size> <out_of_every_n_samples> [window [--benchmark]]\n"
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
"    mono2stereo_s16\n"
"    setbuf <buffer_size>\n"
"    fft_exchange_sides_ff <fft_size>\n"
"    squelch_and_smeter_cc --fifo <squelch_fifo> --outfifo <smeter_fifo> <use_every_nth> <report_every_nth>\n"
"    fifo <buffer_size> <number_of_buffers>\n"
"    invert_u8_u8\n"
"    rtty_line_decoder_u8_u8\n"
"    rtty_baudot2ascii_u8_u8\n"
"    serial_line_decoder_f_u8 <samples_per_bits> [databits [stopbits]]\n"
"    octave_complex_c <samples_to_plot> <out_of_n_samples> [--2d]\n"
"    timing_recovery_cc <algorithm> <decimation> [mu [max_error [--add_q [--output_error | --output_indexes | --octave <show_every_nth> | --octave_save <show_every_nth> <directory> ]]]] \n"
"    psk31_varicode_encoder_u8_u8\n"
"    psk31_varicode_decoder_u8_u8\n"
"    differential_encoder_u8_u8\n"
"    differential_decoder_u8_u8\n"
"    dump_u8\n"
"    psk_modulator_u8_c <n_psk>\n"
"    psk31_interpolate_sine_cc <interpolation>\n"
"    duplicate_samples_ntimes_u8_u8 <sample_size_bytes> <ntimes>\n"
"    bpsk_costas_loop_cc <loop_bandwidth> <damping_factor> [--dd | --decision_directed] [--output_error | --output_dphase | --output_nco | --output_combined <error_file> <dphase_file> <nco_file>]\n"
"    binary_slicer_f_u8\n"
"    simple_agc_cc <rate> [reference [max_gain]]\n"
"    firdes_peak_c <rate> <length> [window [--octave]]\n"
"    peaks_fir_cc <taps_length> [peak_rate × N]\n"
"    repeat_u8 <data_bytes × N>\n"
"    uniform_noise_f\n"
"    gaussian_noise_c\n"
"    awgn_cc <snr_db> [--awgnfile <file>] [--snrshow]\n"
"    pack_bits_8to1_u8_u8\n"
"    pack_bits_1to8_u8_u8\n"
"    firdes_pulse_shaping_filter_f (RRC <samples_per_symbol> <num_taps> <beta> | COSINE <samples_per_symbol>)\n"
"    pulse_shaping_filter_cc (RRC <samples_per_symbol> <num_taps> <beta> | COSINE <samples_per_symbol>)\n"
"    add_n_zero_samples_at_beginning_f <n_zero_samples>\n"
"    generic_slicer_f_u8 <n_symbols>\n"
"    plain_interpolate_cc <interpolation>\n"
"    add_const_cc <i> <q>\n"
"    tee <path> [buffers]\n"
"    pll_cc (1 [alpha] |2 [bandwidth [damping_factor [ko [kd]]]])\n"
"    pattern_search_u8_u8 <values_after> <pattern_values × N>\n" 
"    dbpsk_decoder_c_u8\n" 
"    bfsk_demod_cf <spacing> <filter_length>\n"
"    normalized_timing_variance_u32_f <samples_per_symbol> <initial_sample_offset> [--debug]\n"
"    ?<search_the_function_list>\n"
"    ??<jump_to_function_docs_on_github>\n"
"    =<evaluate_python_expression>\n"
"    shift_addfast_cc <rate>   #only if system supports NEON \n"
"    shift_unroll_cc <rate>\n"
"    logaveragepower_cf <add_db> <fft_size> <avgnumber>\n"
"    fft_one_side_ff <fft_size>\n"
"    convert_f_samplerf <wait_for_this_sample>\n"
"    add_dcoffset_cc\n"
"    fastddc_fwd_cc <decimation> [transition_bw [window]]\n"
"    fastddc_inv_cc <shift_rate> <decimation> [transition_bw [window]]\n"
"    _fft2octave <fft_size>\n"
"    convert_f_i16             #deprecated, use instead: convert_f_s16\n"
"    convert_i16_f             #deprecated, use instead: convert_s16_f\n"
"    floatdump_f               #deprecated, use instead: dump_f\n"
"    mono2stereo_i16           #deprecated, use instead: mono2stereo_s16\n"
"    decode_ima_adpcm_u8_i16   #deprecated, use instead: decode_ima_adpcm_u8_s16\n"
"    encode_ima_adpcm_i16_u8   #deprecated, use instead: encode_ima_adpcm_i16_u8\n"
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
#define TRY_YIELD fflush(stdout);sched_yield()
//unsigned yield_counter=0;

char **argv_global;
int argc_global;

int errhead()
{
    fprintf(stderr, "%s%s%s: ", argv_global[0], ((argc_global>=2)?" ":""), ((argc_global>=2)?argv_global[1]:""));
}

int badsyntax(char* why)
{
    if(why==0) fprintf(stderr, "%s", usage);
    else 
    {
        errhead();
        fprintf(stderr, "%s\n", why);
    }
    return -1;
}

int clipdetect_ff(float* input, int input_size)
{
    for(int i=0;i<input_size;i++)
    {
        if(input[i]<-1.0) { errhead(); fprintf(stderr, "Signal value below -1.0!\n"); return -1; }
        if(input[i]>1.0) { errhead(); fprintf(stderr, "Signal value above 1.0!\n"); return 1; }
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

#define FREAD_U8    fread (input_buffer,    sizeof(unsigned char), the_bufsize, stdin)
#define FWRITE_U8   fwrite (output_buffer,  sizeof(unsigned char), the_bufsize, stdout)
#define FREAD_R     fread (input_buffer,    sizeof(float),      the_bufsize, stdin)
#define FREAD_C     fread (input_buffer,    sizeof(float)*2,    the_bufsize, stdin)
#define FWRITE_R    fwrite (output_buffer,  sizeof(float),      the_bufsize, stdout)
#define FWRITE_C    fwrite (output_buffer,  sizeof(float)*2,    the_bufsize, stdout)
#define FEOF_CHECK  if(feof(stdin)) return 0
//#define BIG_FREAD_C fread(input_buffer, sizeof(float)*2, BIG_BUFSIZE, stdin)
//#define BIG_FWRITE_C fwrite(output_buffer, sizeof(float)*2, BIG_BUFSIZE, stdout)

int init_fifo(int argc, char *argv[])
{
    if(argc>=4)
    {
        if(!strcmp(argv[2],"--fifo"))
        {
            errhead(); fprintf(stderr,"fifo control mode on\n");
            int fd = open(argv[3], O_RDONLY);
            int flags = fcntl(fd, F_GETFL, 0);
            fcntl(fd, F_SETFL, flags | O_NONBLOCK);
            return fd;
        }
        else if(!strcmp(argv[2],"--fd"))  
        {
            //to use this:
            //1. Create a pipe(pipedesc) in your process.
            //2. fork() and execl() your process to run csdr, and give pipedesc[0] as parameter after --fd 
            //  Note: when forking, the child process will get a copy of the file descriptor table! That's why this 
            //  works at all, as file descriptor indexes are normally not transferable between processes, except for a *NIX socket way which is quite complicated... 
            //3. From your parent process, write into pipedesc[1].
            //This is implemented in ddcd, check there to see how to do it!
            int fd;
            if(sscanf(argv[3], "%d",&fd)<=0) return 0;
            errhead();
            fprintf(stderr,"fd control mode on, fd=%d\n", fd);
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



#define UNITROUND_UNIT 4

int unitround(int what)
{
    if(what<=0) return UNITROUND_UNIT;
    return ((what-1)&~(UNITROUND_UNIT-1))+UNITROUND_UNIT;
}

int initialize_buffers()
{
    if(!(the_bufsize=getbufsize())) return 0;
    the_bufsize=unitround(the_bufsize);
    if(env_csdr_print_bufsizes) { errhead(); fprintf(stderr,"buffer size set to %d\n", the_bufsize); }
    input_buffer =  (float*)        malloc(the_bufsize*sizeof(float) * 2); //need the 2× because we might also put complex floats into it
    output_buffer = (float*)        malloc(the_bufsize*sizeof(float) * 2);
    buffer_u8 =     (unsigned char*)malloc(the_bufsize*sizeof(unsigned char));
    buffer_i16 =    (short*)        malloc(the_bufsize*sizeof(short));
    temp_f =        (float*)        malloc(the_bufsize*sizeof(float) * 4);
    if(the_bufsize<=4096) //this is hacky, should be done correctly
    {
        fcntl(STDIN_FILENO, F_SETPIPE_SZ,  4096);
        fcntl(STDOUT_FILENO, F_SETPIPE_SZ, 4096);
    }
    return the_bufsize;
}

int sendbufsize(int size)
{
    if(size<=4096)
    {
        fcntl(STDOUT_FILENO, F_SETPIPE_SZ, 4096);
    }
    //The first word is a preamble, "csdr".
    //If the next csdr process detects it, sets the buffer size according to the second word
    if(!env_csdr_dynamic_bufsize_on) return env_csdr_fixed_bufsize;
    if(env_csdr_print_bufsizes) { errhead(); fprintf(stderr,"next process proposed input buffer size is %d\n", size); }
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
    argc_global=argc;
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

    if(!strcmp(argv[1],"clone") || !strcmp(argv[1],"REM"))
    {
        if(!sendbufsize(initialize_buffers())) return -2;
        clone_(the_bufsize);
    }
#define SET_NONBLOCK(fd) fcntl(fd, F_SETFL, fcntl(fd, F_GETFL, 0) | O_NONBLOCK)

    if(!strcmp(argv[1],"fifo"))
    {
        if(!sendbufsize(initialize_buffers())) return -2;

        int fifo_buffer_size;
        if(argc<=2) return badsyntax("need required parameter (buffer_size)");
        sscanf(argv[2],"%d",&fifo_buffer_size);
        int fifo_num_buffers;
        if(argc<=3) return badsyntax("need required parameter (number of buffers)");
        sscanf(argv[3],"%d",&fifo_num_buffers);

        char** fifo_buffers = (char**)malloc(sizeof(char*)*fifo_num_buffers);
        for(int i=0;i<fifo_num_buffers;i++) fifo_buffers[i]=(char*)malloc(sizeof(char)*fifo_buffer_size);

        SET_NONBLOCK(STDIN_FILENO);
        SET_NONBLOCK(STDOUT_FILENO);

        fd_set read_fds;
        FD_ZERO(&read_fds);
        FD_SET(STDIN_FILENO, &read_fds);
        fd_set write_fds;
        FD_ZERO(&write_fds);
        FD_SET(STDOUT_FILENO, &write_fds);

        int highfd = ((STDOUT_FILENO > STDIN_FILENO) ? STDOUT_FILENO : STDIN_FILENO) + 1;

        int fifo_actual_buffer_wr = fifo_num_buffers - 1;
        int fifo_actual_buffer_rd = 0;
        int fifo_actual_buffer_wr_pos = 0;
        int fifo_actual_buffer_rd_pos = 0;
        int fifo_error = 0;
        int fifo_overrun_shown = 0;

        for(;;)
        {
            select(highfd, &read_fds, NULL, NULL, NULL);

            //try to read until buffer is full
            if(FD_ISSET(STDIN_FILENO, &read_fds)) for(;;)
            {
                int read_bytes=read(STDIN_FILENO, fifo_buffers[fifo_actual_buffer_rd]+fifo_actual_buffer_rd_pos, fifo_buffer_size-fifo_actual_buffer_rd_pos);
                //fprintf(stderr, "r %d %d | %d %d\n", read_bytes, fifo_buffer_size-fifo_actual_buffer_rd_pos, fifo_actual_buffer_rd, fifo_actual_buffer_rd_pos);
                if(!read_bytes || ((read_bytes<0)&&(fifo_error=read_bytes)) ) break;
                fifo_actual_buffer_rd_pos+=read_bytes;
                if(!((fifo_actual_buffer_rd==fifo_actual_buffer_wr-1)||(fifo_actual_buffer_wr==0&&fifo_actual_buffer_rd==fifo_num_buffers-1)))
                {
                    if(fifo_actual_buffer_rd_pos==fifo_buffer_size)
                    {
                        fifo_overrun_shown = 0;
                        fifo_actual_buffer_rd++;
                        fifo_actual_buffer_rd_pos = 0;
                        if(fifo_actual_buffer_rd>=fifo_num_buffers) fifo_actual_buffer_rd=0;
                    }
                }
                else
                {
                    if(fifo_actual_buffer_rd_pos==fifo_buffer_size)
                    {
                        fifo_actual_buffer_rd_pos = 0; //rewrite same buffer
                        if(!fifo_overrun_shown) { fifo_overrun_shown=1; errhead(); fprintf(stderr, "circular buffer full, dropping samples\n"); }
                    }
                }
            }
            //try to write until buffer is empty
            if(FD_ISSET(STDOUT_FILENO, &write_fds)) for(;;)
            {
                if(fifo_actual_buffer_wr == fifo_actual_buffer_rd) break;
                int written_bytes=write(STDOUT_FILENO, fifo_buffers[fifo_actual_buffer_wr]+fifo_actual_buffer_wr_pos, fifo_buffer_size-fifo_actual_buffer_wr_pos);
                //fprintf(stderr, "w %d %d | %d %d\n", written_bytes, fifo_buffer_size-fifo_actual_buffer_wr_pos, fifo_actual_buffer_wr, fifo_actual_buffer_wr_pos);
                if(!written_bytes || ((written_bytes<0)&&(fifo_error=written_bytes)) ) break;
                fifo_actual_buffer_wr_pos+=written_bytes;
                if(fifo_actual_buffer_wr_pos==fifo_buffer_size)
                {
                    fifo_actual_buffer_wr++;
                    fifo_actual_buffer_wr_pos = 0;
                    if(fifo_actual_buffer_wr>=fifo_num_buffers) fifo_actual_buffer_wr=0;
                }

            }
            if(fifo_error&&errno!=11) { errhead(); fprintf(stderr,"fifo_error (%d)", errno); return -1; }
        }

        return -1;

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
    if((!strcmp(argv[1],"convert_i16_f")) || (!strcmp(argv[1],"convert_s16_f")))
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
    if(!strcmp(argv[1],"convert_f_s24"))
    {
        int bigendian = (argc>2) && (!strcmp(argv[2],"--bigendian"));
        unsigned char* s24buffer = (unsigned char*)malloc(sizeof(unsigned char)*the_bufsize*3);
        if(!sendbufsize(initialize_buffers())) return -2;
        for(;;)
        {
            FEOF_CHECK;
            FREAD_R;
            convert_f_s24(input_buffer, s24buffer, the_bufsize, bigendian);
            fwrite(s24buffer, sizeof(unsigned char)*3, the_bufsize, stdout);
            TRY_YIELD;
        }
    }
    if(!strcmp(argv[1],"convert_s24_f"))
    {
        int bigendian = (argc>2) && (!strcmp(argv[2],"--bigendian"));
        unsigned char* s24buffer = (unsigned char*)malloc(sizeof(unsigned char)*the_bufsize*3);
        if(!sendbufsize(initialize_buffers())) return -2;
        for(;;)
        {
            FEOF_CHECK;
            fread(s24buffer, sizeof(unsigned char)*3, the_bufsize, stdin);
            convert_s24_f(s24buffer, output_buffer, the_bufsize, bigendian);
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
            starting_phase=shift_math_cc((complexf*)input_buffer, (complexf*)output_buffer, the_bufsize, rate, starting_phase);
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
        errhead();
        fprintf(stderr,"LUT initialized\n");
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

    if(!strcmp(argv[1],"shift_addfast_cc"))
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
            shift_addfast_data_t data=shift_addfast_init(rate);
            errhead();
            fprintf(stderr,"reinitialized to %g\n",rate);
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
                    starting_phase=shift_addfast_cc((complexf*)ibufptr, (complexf*)obufptr, current_size, &data, starting_phase);
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


    if(!strcmp(argv[1],"shift_unroll_cc"))
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
            shift_unroll_data_t data=shift_unroll_init(rate, 1024);
            errhead();
            fprintf(stderr,"reinitialized to %g\n",rate);
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
                    starting_phase=shift_unroll_cc((complexf*)ibufptr, (complexf*)obufptr, current_size, &data, starting_phase);
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
            errhead();
            fprintf(stderr,"reinitialized to %g\n",rate);
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
        errhead(); fprintf(stderr,"tau = %g, sample_rate = %d\n",tau,sample_rate);
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
            if(nan_detect) { errhead(); fprintf(stderr, "NaN detected!\n"); }
            fwrite(input_buffer, sizeof(float), the_bufsize, stdout);
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"floatdump_f") || !strcmp(argv[1],"dump_f"))
    {
        if(!sendbufsize(initialize_buffers())) return -2;
        for(;;)
        {
            FEOF_CHECK;
            FREAD_R;
            for(int i=0; i<the_bufsize;i++) printf("%g ",input_buffer[i]);
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

        int taps_length=firdes_filter_len(transition_bw);
        fprintf(stderr,"fir_decimate_cc: taps_length = %d\n",taps_length);

        while (env_csdr_fixed_big_bufsize < taps_length*2) env_csdr_fixed_big_bufsize*=2; //temporary fix for buffer size if [transition_bw] is low
        //fprintf(stderr, "env_csdr_fixed_big_bufsize = %d\n", env_csdr_fixed_big_bufsize);

        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize/factor);


        int padded_taps_length = taps_length;
        float *taps;
#define NEON_ALIGNMENT (4*4*2)
#ifdef NEON_OPTS
        errhead(); fprintf(stderr,"taps_length = %d\n", taps_length);
        padded_taps_length = taps_length+(NEON_ALIGNMENT/4)-1 - ((taps_length+(NEON_ALIGNMENT/4)-1)%(NEON_ALIGNMENT/4));
        errhead(); fprintf(stderr,"padded_taps_length = %d\n", padded_taps_length);

        taps = (float*) (float*)malloc((padded_taps_length+NEON_ALIGNMENT)*sizeof(float));
        errhead(); fprintf(stderr,"taps = %x\n", taps);
        taps =  (float*)((((unsigned)taps)+NEON_ALIGNMENT-1) & ~(NEON_ALIGNMENT-1));
        errhead(); fprintf(stderr,"NEON aligned taps = %x\n", taps);
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
            TRY_YIELD;
            input_skip=factor*output_size;
            memmove((complexf*)input_buffer,((complexf*)input_buffer)+input_skip,(the_bufsize-input_skip)*sizeof(complexf)); //memmove lets the source and destination overlap
            fread(((complexf*)input_buffer)+(the_bufsize-input_skip), sizeof(complexf), input_skip, stdin);
            //fprintf(stderr,"iskip=%d output_size=%d start=%x target=%x skipcount=%x \n",input_skip,output_size,input_buffer, ((complexf*)input_buffer)+(BIG_BUFSIZE-input_skip),(BIG_BUFSIZE-input_skip));
        }
    }

    if(!strcmp(argv[1],"fir_interpolate_cc"))
    {
        bigbufs=1;

        if(argc<=2) return badsyntax("need required parameter (interpolation factor)");

        int factor;
        sscanf(argv[2],"%d",&factor);
        assert(factor >= 1);

        float transition_bw = 0.05;
        if(argc>=4) sscanf(argv[3],"%g",&transition_bw);
        assert(transition_bw >= 0 && transition_bw < 1.);

        window_t window = WINDOW_DEFAULT;
        if(argc>=5)
        {
            window=firdes_get_window_from_string(argv[4]);
        }
        else {errhead(); fprintf(stderr,"window = %s\n",firdes_get_string_from_window(window));}

        int taps_length=firdes_filter_len(transition_bw);
        errhead(); fprintf(stderr,"taps_length = %d\n",taps_length);
        assert(taps_length > 0);

        while (env_csdr_fixed_big_bufsize < taps_length*2) env_csdr_fixed_big_bufsize*=2; //temporary fix for buffer size if [transition_bw] is low
        //fprintf(stderr, "env_csdr_fixed_big_bufsize = %d\n", env_csdr_fixed_big_bufsize);

        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize*factor);
        assert(the_bufsize > 0);

        float *taps;
        taps=(float*)malloc(taps_length*sizeof(float));
        assert(taps);

        firdes_lowpass_f(taps,taps_length,0.5/(float)factor,window);

        int input_skip=0;
        int output_size=0;
        float* interp_output_buffer = (float*)malloc(sizeof(float)*2*the_bufsize*factor);
        for(;;)
        {
            FEOF_CHECK;
            output_size=fir_interpolate_cc((complexf*)input_buffer, (complexf*)interp_output_buffer, the_bufsize, factor, taps, taps_length);
            //fprintf(stderr, "os %d\n",output_size);
            fwrite(interp_output_buffer, sizeof(complexf), output_size, stdout);
            TRY_YIELD;
            input_skip=output_size/factor;
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
        else { errhead(); fprintf(stderr,"window = %s\n",firdes_get_string_from_window(window)); }

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
        else { errhead(); fprintf(stderr,"window = %s\n",firdes_get_string_from_window(window));} 

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
        else { errhead(); fprintf(stderr,"window = %s\n",firdes_get_string_from_window(window)); }

        if(suboptimal) { errhead(); fprintf(stderr,"note: suboptimal rational resampler chosen.\n"); }

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

        int num_poly_points = 12;
        if(argc>=4) sscanf(argv[3],"%d",&num_poly_points);
        if(num_poly_points&1) return badsyntax("num_poly_points should be even");
        if(num_poly_points<2) return badsyntax("num_poly_points should be >= 2");

        int use_prefilter = 0;
        float transition_bw=0.03;
        window_t window = WINDOW_DEFAULT;
        if(argc>=5)
        {
            if(!strcmp(argv[4], "--prefilter")) 
            {
                errhead(); fprintf(stderr, "using prefilter with default values\n"); 
                use_prefilter = 1;
            }
            else 
            {
                sscanf(argv[4],"%g",&transition_bw);
                if(argc>=6) window = firdes_get_window_from_string(argv[5]);
            }
        }
        errhead(); fprintf(stderr,"use_prefilter = %d, num_poly_points = %d, transition_bw = %g, window = %s\n", 
            use_prefilter, num_poly_points, transition_bw, firdes_get_string_from_window(window));

        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize / rate);

        if(rate==1) clone_(the_bufsize); //copy input to output in this special case (and stick in this function).

        //Generate filter taps
        int taps_length = 0;
        float* taps = NULL;
        if(use_prefilter)
        {
            taps_length = firdes_filter_len(transition_bw);
            errhead(); fprintf(stderr,"taps_length = %d\n",taps_length);
            taps = (float*)malloc(sizeof(float)*taps_length);
            firdes_lowpass_f(taps, taps_length, 0.5/(rate-transition_bw), window); //0.6 const to compensate rolloff
            //for(int=0;i<taps_length; i++) fprintf(stderr,"%g ",taps[i]);
        }
        else { errhead(); fprintf(stderr,"not using taps\n"); }
        fractional_decimator_ff_t d = fractional_decimator_ff_init(rate, num_poly_points, taps, taps_length); 
        for(;;)
        {
            FEOF_CHECK;
            if(d.input_processed==0) d.input_processed=the_bufsize;
            else memcpy(input_buffer, input_buffer+d.input_processed, sizeof(float)*(the_bufsize-d.input_processed));
            fread(input_buffer+(the_bufsize-d.input_processed), sizeof(float), d.input_processed, stdin);
            fractional_decimator_ff(input_buffer, output_buffer, the_bufsize, &d);
            fwrite(output_buffer, sizeof(float), d.output_size, stdout);
            //fprintf(stderr, "os = %d, ip = %d\n", d.output_size, d.input_processed);
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"old_fractional_decimator_ff"))
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
        else { errhead(); fprintf(stderr,"window = %s\n",firdes_get_string_from_window(window)); }

        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize / rate);

        if(rate==1) clone_(the_bufsize); //copy input to output in this special case (and stick in this function).

        //Generate filter taps
        int taps_length = firdes_filter_len(transition_bw);
        errhead(); fprintf(stderr,"taps_length = %d\n",taps_length); 
        float* taps = (float*)malloc(sizeof(float)*taps_length);
        firdes_lowpass_f(taps, taps_length, 0.59*0.5/(rate-transition_bw), window); //0.6 const to compensate rolloff
        //for(int=0;i<taps_length; i++) fprintf(stderr,"%g ",taps[i]);

        static old_fractional_decimator_ff_t d; //in .bss => initialized to zero
        for(;;)
        {
            FEOF_CHECK;
            if(d.input_processed==0) d.input_processed=the_bufsize;
            else memcpy(input_buffer, input_buffer+d.input_processed, sizeof(float)*(the_bufsize-d.input_processed));
            fread(input_buffer+(the_bufsize-d.input_processed), sizeof(float), d.input_processed, stdin);
            d = old_fractional_decimator_ff(input_buffer, output_buffer, the_bufsize, rate, taps, taps_length, d);
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
        if(benchmark) { errhead(); fprintf(stderr,"benchmarking..."); }
        FFT_PLAN_T* plan=make_fft_c2c(fft_size, windowed, output, 1, benchmark);
        if(benchmark) fprintf(stderr," done\n");
        if(octave) printf("setenv(\"GNUTERM\",\"X11 noraise\");y=zeros(1,%d);semilogy(y,\"ydatasource\",\"y\");\n",fft_size);
        float *windowt;
        windowt = precalculate_window(fft_size, window);
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
            //apply_window_c(input,windowed,fft_size,window);
            apply_precalculated_window_c(input,windowed,fft_size,windowt);
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

    if(!strcmp(argv[1],"logaveragepower_cf"))
    {
        bigbufs=1;
        if(argc<=4) return badsyntax("need required parameters (add_db, fft_size, avgnumber)"); 
        float add_db=0;
        int avgnumber=0;
        int fft_size=0;
        
        sscanf(argv[2],"%g",&add_db);
        sscanf(argv[3],"%d",&fft_size);
        sscanf(argv[4],"%d",&avgnumber);
        
        float *input = malloc(sizeof(float)*2 * fft_size);
        float *output = malloc(sizeof(float) * fft_size);

        add_db -= 10.0*log10(avgnumber);
        for(;;)
        {
            int i,n;
            for(i = 0; i < fft_size; i++) {
                output[i] = 0;
            }
            FEOF_CHECK;
            for(n = 0; n < avgnumber; n++) {
                fread (input, sizeof(float)*2, fft_size, stdin);
                accumulate_power_cf((complexf*)input, output, fft_size);
            }
            log_ff(output, output, fft_size, add_db);
            fwrite (output, sizeof(float), fft_size, stdout);
            TRY_YIELD;
        }
        return 0;
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

    if(!strcmp(argv[1],"fft_one_side_ff"))
    {
        if(argc<=2) return badsyntax("need required parameters (fft_size)");
        int fft_size;
        sscanf(argv[2],"%d",&fft_size);
        if(!getbufsize()) return -2; 
        sendbufsize(fft_size);
        float* input_buffer_s1 = (float*)malloc(sizeof(float)*fft_size/2);
        float* input_buffer_s2 = (float*)malloc(sizeof(float)*fft_size/2);
        for(;;)
        {
            FEOF_CHECK;
            fread(input_buffer_s1, sizeof(float), fft_size/2, stdin);
            fread(input_buffer_s2, sizeof(float), fft_size/2, stdin);
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

    if(!strcmp(argv[1],"fft_benchmark"))
    {
        if(argc<=3) return badsyntax("need required parameters (fft_size, fft_cycles)");
        int fft_size;
        sscanf(argv[2],"%d",&fft_size);
        int fft_cycles;
        sscanf(argv[3],"%d",&fft_cycles);

        int benchmark=(argc>=5)&&!strcmp(argv[4],"--benchmark");
        errhead(); fprintf(stderr,"FFT library used: %s\n",FFT_LIBRARY_USED);

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
        errhead(); fprintf(stderr,"initializing... ");
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
        errhead(); fprintf(stderr,"%d transforms of %d processed in %g seconds, %g seconds each.\n",fft_cycles,fft_size,time_taken_fft,time_taken_fft/fft_cycles);
        return 0;
    }

    if(!strcmp(argv[1],"bandpass_fir_fft_cc")) //this command does not exist as a separate function
    {
        float low_cut;
        float high_cut;
        float transition_bw;
        window_t window = WINDOW_DEFAULT;
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
        if(argc>=6) window=firdes_get_window_from_string(argv[5]);
        else { errhead(); fprintf(stderr,"window = %s\n",firdes_get_string_from_window(window)); }

        //calculate the FFT size and the other length parameters
        int taps_length=firdes_filter_len(transition_bw); //the number of non-zero taps
        int fft_size=next_pow2(taps_length); //we will have to pad the taps with zeros until the next power of 2 for FFT
        //the number of padding zeros is the number of output samples we will be able to take away after every processing step, and it looks sane to check if it is large enough.
        if (fft_size-taps_length<200) fft_size<<=1;
        int input_size = fft_size - taps_length + 1;
        int overlap_length = taps_length - 1;
        errhead(); fprintf(stderr,"(fft_size = %d) = (taps_length = %d) + (input_size = %d) - 1\n(overlap_length = %d) = taps_length - 1\n", fft_size, taps_length, input_size, overlap_length );
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
            errhead(); fprintf(stderr,"filter initialized, low_cut = %g, high_cut = %g\n",low_cut,high_cut);
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
        errhead(); fprintf(stderr, "flowcontrol_bufsize = %d, flowcontrol_sleep = %d\n", flowcontrol_bufsize, flowcontrol_sleep);
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
                    fprintf( stderr, "through: %lu bytes/s, buffer #%d\n", (unsigned long)floor((float)buffer_count*the_bufsize*sizeof(float)/timetaken), buffer_count );
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
        errhead(); fprintf(stderr, "initial squelch level is %g\n", squelch_level);
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
            if(read_fifo_ctl(fd,"%g\n",&squelch_level)) { errhead(); fprintf(stderr, "new squelch level is %g\n", squelch_level); }
            TRY_YIELD;
        }
    }

    /*
      ______        _   _____  _____   _____ 
     |  ____|      | | |  __ \|  __ \ / ____|
     | |__ __ _ ___| |_| |  | | |  | | |     
     |  __/ _` / __| __| |  | | |  | | |     
     | | | (_| \__ \ |_| |__| | |__| | |____ 
     |_|  \__,_|___/\__|_____/|_____/ \_____|

    */                                         
                                         
    if( !strcmp(argv[1],"fastddc_fwd_cc") ) //<decimation> [transition_bw [window]]
    {   
        
        int decimation;
        if(argc<=2) return badsyntax("need required parameter (decimation)");
        sscanf(argv[2],"%d",&decimation);
        
        float transition_bw = 0.05;
        if(argc>3) sscanf(argv[3],"%g",&transition_bw);

        window_t window = WINDOW_DEFAULT;
        if(argc>4)  window=firdes_get_window_from_string(argv[4]);
        else { errhead(); fprintf(stderr,"window = %s\n",firdes_get_string_from_window(window)); }

        fastddc_t ddc; 
        if(fastddc_init(&ddc, transition_bw, decimation, 0)) { badsyntax("error in fastddc_init()"); return 1; }
        fastddc_print(&ddc,"fastddc_fwd_cc");

        if(!initialize_buffers()) return -2;
        sendbufsize(ddc.fft_size);

        //make FFT plan
        complexf* input =    (complexf*)fft_malloc(sizeof(complexf)*ddc.fft_size);
        complexf* windowed = (complexf*)fft_malloc(sizeof(complexf)*ddc.fft_size);
        complexf* output =   (complexf*)fft_malloc(sizeof(complexf)*ddc.fft_size);

        for(int i=0;i<ddc.fft_size;i++) iof(input,i)=qof(input,i)=0; //null the input buffer

        int benchmark = 1; 
        if(benchmark) { errhead(); fprintf(stderr,"benchmarking FFT..."); }
        FFT_PLAN_T* plan=make_fft_c2c(ddc.fft_size, windowed, output, 1, benchmark);
        if(benchmark) fprintf(stderr," done\n");

        for(;;)
        {
            FEOF_CHECK;
            //overlapped FFT
            for(int i=0;i<ddc.overlap_length;i++) input[i]=input[i+ddc.input_size];
            fread(input+ddc.overlap_length, sizeof(complexf), ddc.input_size, stdin);
            //apply_window_c(input,windowed,ddc.fft_size,window);
            memcpy(windowed, input, ddc.fft_size*sizeof(complexf)); //we can switch off windows; TODO: it is likely that we shouldn't apply a window to both the FFT and the filter.
            fft_execute(plan);
            fwrite(output, sizeof(complexf), ddc.fft_size, stdout);
            TRY_YIELD;
        }
    }

    if( !strcmp(argv[1],"fastddc_inv_cc") ) //<shift_rate> <decimation> [transition_bw [window]]
    {   
        float shift_rate;
        int plusarg=0;

        int fd;
        if(fd=init_fifo(argc,argv))
        {
            while(!read_fifo_ctl(fd,"%g\n",&shift_rate)) usleep(10000);
            plusarg=1;
        }
        else
        {
            if(argc<=2) return badsyntax("need required parameter (rate)"); 
            sscanf(argv[2],"%g",&shift_rate);
        }

        int decimation;
        if(argc<=3+plusarg) return badsyntax("need required parameter (decimation)");
        sscanf(argv[3+plusarg],"%d",&decimation);
        //fprintf(stderr, "dec=%d %d\n", decimation);

        float transition_bw = 0.05;
        if(argc>4+plusarg) sscanf(argv[4+plusarg],"%g",&transition_bw);

        window_t window = WINDOW_DEFAULT;
        if(argc>5+plusarg)  window=firdes_get_window_from_string(argv[5+plusarg]);
        else { errhead(); fprintf(stderr,"window = %s\n",firdes_get_string_from_window(window)); }

        for(;;)
        {

        fastddc_t ddc; 
        if(fastddc_init(&ddc, transition_bw, decimation, shift_rate)) { badsyntax("error in fastddc_init()"); return 1; }
        fastddc_print(&ddc,"fastddc_inv_cc");

        if(!initialize_buffers()) return -2;
        sendbufsize(ddc.post_input_size/ddc.post_decimation); //TODO not exactly correct

        //prepare making the filter and doing FFT on it
        complexf* taps=(complexf*)calloc(sizeof(complexf),ddc.fft_size); //initialize to zero
        complexf* taps_fft=(complexf*)malloc(sizeof(complexf)*ddc.fft_size);
        FFT_PLAN_T* plan_taps = make_fft_c2c(ddc.fft_size, taps, taps_fft, 1, 0); //forward, don't benchmark (we need this only once)

        //make the filter
        float filter_half_bw = 0.5/decimation;
        errhead(); fprintf(stderr, "preparing a bandpass filter of [%g, %g] cutoff rates. Real transition bandwidth is: %g\n", (-shift_rate)-filter_half_bw, (-shift_rate)+filter_half_bw, 4.0/ddc.taps_length);
        firdes_bandpass_c(taps, ddc.taps_length, (-shift_rate)-filter_half_bw, (-shift_rate)+filter_half_bw, window);
        fft_execute(plan_taps);
        fft_swap_sides(taps_fft,ddc.fft_size);

        //make FFT plan
        complexf* inv_input =    (complexf*)fft_malloc(sizeof(complexf)*ddc.fft_inv_size);
        complexf* inv_output =   (complexf*)fft_malloc(sizeof(complexf)*ddc.fft_inv_size);
        errhead(); fprintf(stderr,"benchmarking FFT...");
        FFT_PLAN_T* plan_inverse = make_fft_c2c(ddc.fft_inv_size, inv_input, inv_output, 0, 1); //inverse, do benchmark
        fprintf(stderr," done\n");
        
        //alloc. buffers
        complexf* input =    (complexf*)fft_malloc(sizeof(complexf)*ddc.fft_size);
        complexf* output =   (complexf*)fft_malloc(sizeof(complexf)*ddc.post_input_size);

        decimating_shift_addition_status_t shift_stat;
        bzero(&shift_stat, sizeof(shift_stat));
        for(;;)
        {
            FEOF_CHECK;
            fread(input, sizeof(complexf), ddc.fft_size, stdin);
            shift_stat = fastddc_inv_cc(input, output, &ddc, plan_inverse, taps_fft, shift_stat);
            fwrite(output, sizeof(complexf), shift_stat.output_size, stdout);
            //fprintf(stderr, "ss os = %d\n", shift_stat.output_size);
            TRY_YIELD;
            if(read_fifo_ctl(fd,"%g\n",&shift_rate)) break;
        }

        }
    }

    if( !strcmp(argv[1], "_fft2octave") ) 
    {
        int fft_size;
        if(argc<=2) return badsyntax("need required parameter (fft_size)");
        sscanf(argv[2],"%d",&fft_size);

        complexf* fft_input=(complexf*)malloc(sizeof(complexf)*fft_size);
        initialize_buffers();
        if(!sendbufsize(fft_size)) return -2;

        printf("setenv(\"GNUTERM\",\"X11 noraise\");y=zeros(1,%d);semilogy(y,\"ydatasource\",\"y\");\n",fft_size);
        for(;;)
        {
            FEOF_CHECK;
            fread(fft_input, sizeof(complexf), fft_size, stdin);
            printf("fftdata=[");
            //we have to swap the two parts of the array to get a valid spectrum
            for(int i=fft_size/2;i<fft_size;i++) printf("(%g)+(%g)*i ",iof(fft_input,i),qof(fft_input,i));
            for(int i=0;i<fft_size/2;i++) printf("(%g)+(%g)*i ",iof(fft_input,i),qof(fft_input,i)); 
            printf(
                "];\n"
                "y=abs(fftdata);\n"
                "refreshdata;\n"
            );
        }
    }

    /*
      _____  _       _ _        _                       _           
     |  __ \(_)     (_) |      | |                     | |          
     | |  | |_  __ _ _| |_ __ _| |  _ __ ___   ___   __| | ___  ___ 
     | |  | | |/ _` | | __/ _` | | | '_ ` _ \ / _ \ / _` |/ _ \/ __|
     | |__| | | (_| | | || (_| | | | | | | | | (_) | (_| |  __/\__ \
     |_____/|_|\__, |_|\__\__,_|_| |_| |_| |_|\___/ \__,_|\___||___/
                __/ |                                               
               |___/                                                
    */

    if(!strcmp(argv[1],"psk31_varicode_decoder_u8_u8"))
    {
        unsigned long long status_shr = 0;
        unsigned char output;
        if(!sendbufsize(initialize_buffers())) return -2;
        unsigned char i=0;
        for(;;)
        {
            if((output=psk31_varicode_decoder_push(&status_shr, getchar()))) { putchar(output); fflush(stdout); }
            if(i++) continue; //do the following at every 256th execution of the loop body:
            FEOF_CHECK;
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"invert_u8_u8"))
    {
        if(!sendbufsize(initialize_buffers())) return -2;
        unsigned char i=0;
        for(;;)
        {
            putchar(!getchar());
            if(i++) continue; //do the following at every 256th execution of the loop body:
            FEOF_CHECK;
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"rtty_line_decoder_u8_u8"))
    {
        static rtty_baudot_decoder_t status_baudot; //created on .bss -> initialized to 0
        unsigned char output;
        if(!sendbufsize(initialize_buffers())) return -2;
        unsigned char i=0;
        for(;;)
        {
            if((output=rtty_baudot_decoder_push(&status_baudot, getchar()))) { putchar(output); fflush(stdout); }
            if(i++) continue; //do the following at every 256th execution of the loop body:
            FEOF_CHECK;
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"rtty_baudot2ascii_u8_u8"))
    {
        unsigned char fig_mode = 0;
        unsigned char output;
        if(!sendbufsize(initialize_buffers())) return -2;
        unsigned char i=0;
        for(;;)
        {
            if((output=rtty_baudot_decoder_lookup(&fig_mode, getchar()))) { putchar(output); fflush(stdout); }
            if(i++) continue; //do the following at every 256th execution of the loop body:
            FEOF_CHECK;
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"binary_slicer_f_u8"))
    {
        if(!sendbufsize(initialize_buffers())) return -2;
        for(;;)
        {
            FEOF_CHECK;
            if(!FREAD_R) break;
            binary_slicer_f_u8(input_buffer, (unsigned char*)output_buffer, the_bufsize);
            FWRITE_U8;
            TRY_YIELD;
        }
        return 0;
    }

    if(!strcmp(argv[1],"serial_line_decoder_f_u8")) //<samples_per_bits> [databits [stopbits]]
    {
        bigbufs=1;

        serial_line_t serial;

        if(argc<=2) return badsyntax("need required parameter (samples_per_bits)");
        sscanf(argv[2],"%f",&serial.samples_per_bits);
        if(serial.samples_per_bits<1) return badsyntax("samples_per_bits should be at least 1.");
        if(serial.samples_per_bits<5) fprintf(stderr, "%s: warning: this algorithm does not work well if samples_per_bits is too low. It should be at least 5.\n", argv[1]);

        serial.databits=8;
        if(argc>3) sscanf(argv[3],"%d",&serial.databits);
        if(serial.databits>8 || serial.databits<1) return badsyntax("databits should be between 1 and 8.");

        serial.stopbits=1;
        if(argc>4) sscanf(argv[4],"%f",&serial.stopbits);
        if(serial.stopbits<1) return badsyntax("stopbits should be equal or above 1.");

        serial.bit_sampling_width_ratio = 0.4;
        serial.input_used=0;

        if(!sendbufsize(initialize_buffers())) return -2;

        for(;;)
        {
            FEOF_CHECK;
            if(serial.input_used)
            {
                memmove(input_buffer, input_buffer+serial.input_used, sizeof(float)*(the_bufsize-serial.input_used));
                fread(input_buffer+(the_bufsize-serial.input_used), sizeof(float), serial.input_used, stdin);
            }
            else fread(input_buffer, sizeof(float), the_bufsize, stdin); //should happen only on the first run
            serial_line_decoder_f_u8(&serial,input_buffer, (unsigned char*)output_buffer, the_bufsize);
            //printf("now in | ");
            if(serial.input_used==0) { errhead(); fprintf(stderr, "error: serial_line_decoder_f_u8() got stuck.\n"); return -3; }
            //printf("now out %d | ", serial.output_size);
            fwrite(output_buffer, sizeof(unsigned char), serial.output_size, stdout);
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"pll_cc"))
    {
        pll_t pll;

        if(argc<=2) return badsyntax("need required parameter (pll_type)");
        sscanf(argv[2],"%d",(int*)&pll.pll_type);
        if(pll.pll_type == PLL_P_CONTROLLER)
        {
                float alpha = 0.01;
                if(argc>3) sscanf(argv[3],"%f",&alpha);
                pll_cc_init_p_controller(&pll, alpha);
        }
        else if(pll.pll_type == PLL_PI_CONTROLLER)
        {
            float bandwidth = 0.01, ko = 10, kd=0.1, damping_factor = 0.707;
            if(argc>3) sscanf(argv[3],"%f",&bandwidth);
            if(argc>4) sscanf(argv[4],"%f",&damping_factor);
            if(argc>5) sscanf(argv[5],"%f",&ko);
            if(argc>6) sscanf(argv[6],"%f",&kd);
            pll_cc_init_pi_controller(&pll, bandwidth, ko, kd, damping_factor);
            errhead(); fprintf(stderr, "bw=%f damping=%f ko=%f kd=%f alpha=%f beta=%f\n", bandwidth, damping_factor, ko, kd, pll.alpha, pll.beta);
            //  pll.filter_taps_a[0], pll.filter_taps_a[1], pll.filter_taps_a[2], pll.filter_taps_b[0], pll.filter_taps_b[1], pll.filter_taps_b[2]);
        }
        else return badsyntax("invalid pll_type. Valid values are:\n\t1: PLL_P_CONTROLLER\n\t2: PLL_PI_CONTROLLER");

        if(!sendbufsize(initialize_buffers())) return -2;

        for(;;)
        {
            FEOF_CHECK;
            FREAD_C;
            //fprintf(stderr, "| i");
            // pll_cc(&pll, (complexf*)input_buffer, output_buffer, NULL, the_bufsize);
            // fwrite(output_buffer, sizeof(float), the_bufsize, stdout);
            pll_cc(&pll, (complexf*)input_buffer, NULL, (complexf*)output_buffer, the_bufsize);
            fwrite(output_buffer, sizeof(complexf), the_bufsize, stdout);
            //fprintf(stderr, "| o");
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"timing_recovery_cc")) //<algorithm> <decimation> [mu [max_error [--add_q [--output_error | --output_indexes | --octave <show_every_nth> | --octave_save <show_every_nth> <directory> ]]]] \n"
    {
        if(argc<=2) return badsyntax("need required parameter (algorithm)");
        timing_recovery_algorithm_t algorithm = timing_recovery_get_algorithm_from_string(argv[2]);
        //if(algorithm == TIMING_RECOVERY_ALGORITHM_DEFAULT) 
        //  fprintf(stderr,"#timing_recovery_cc: algorithm = %s\n",timing_recovery_get_string_from_algorithm(algorithm));
        if(argc<=3) return badsyntax("need required parameter (decimation factor)");
        int decimation;
        sscanf(argv[3],"%d",&decimation);
        if(decimation<=4 || decimation&3) return badsyntax("decimation factor should be a positive integer divisible by 4");

        float loop_gain = 0.5;
        if(argc>4) sscanf(argv[4],"%f",&loop_gain);

        float max_error = 2;
        if(argc>5) sscanf(argv[5],"%f",&max_error);

        int add_q = !!(argc>=7 && !strcmp(argv[6], "--add_q"));

        int debug_every_nth = -1;
        int output_error = 0;
        int output_indexes = 0;
        int octave_save = 0;
        char* octave_save_path = NULL;
        if(argc>=8+add_q && (!strcmp(argv[6+add_q], "--octave") || (octave_save = !strcmp(argv[6+add_q], "--octave_save")))) 
        {
            debug_every_nth = atoi(argv[7+add_q]);
            if(debug_every_nth<0) return badsyntax("debug_every_nth should be >= 0");
        }
        if(octave_save)
        {
            if(argc>=9+add_q) octave_save_path = argv[8+add_q]; 
            else octave_save_path = "figs";
        }
        if(debug_every_nth<0) { errhead(); fprintf(stderr, "--add_q mode on\n"); }

        if(argc>=(7+add_q) && !strcmp(argv[6+add_q], "--output_error")) output_error = 1;
        float* timing_error = NULL;
        if(output_error) timing_error = (float*)malloc(sizeof(float)*the_bufsize);
        if(output_error) { errhead(); fprintf(stderr, "--output_error mode\n"); }

        if(argc>=(7+add_q) && !strcmp(argv[6+add_q], "--output_indexes")) output_indexes = 1;
        unsigned* sampled_indexes = NULL;
        if(output_indexes) sampled_indexes = (unsigned*)malloc(sizeof(float)*the_bufsize);
        if(output_indexes) { errhead(); fprintf(stderr, "--output_indexes mode\n"); }

        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize/decimation);

        timing_recovery_state_t state = timing_recovery_init(algorithm, decimation, add_q, loop_gain, max_error, debug_every_nth, octave_save_path);

        FREAD_C;
        unsigned buffer_start_counter = 0;
        for(;;)
        {
            FEOF_CHECK;
            timing_recovery_cc((complexf*)input_buffer, (complexf*)output_buffer, the_bufsize, timing_error, (int*)sampled_indexes, &state);
            //fprintf(stderr, "trcc is=%d, os=%d, ip=%d\n",the_bufsize, state.output_size, state.input_processed);
            if(timing_error) fwrite(timing_error, sizeof(float), state.output_size, stdout);
            else if(sampled_indexes) 
            {
                for(int i=0;i<state.output_size;i++) sampled_indexes[i]+=buffer_start_counter;
                fwrite(sampled_indexes, sizeof(unsigned), state.output_size, stdout);
            }
            else fwrite(output_buffer, sizeof(complexf), state.output_size, stdout);
            TRY_YIELD;
            //fprintf(stderr, "state.input_processed = %d\n", state.input_processed);
            buffer_start_counter+=state.input_processed;
            memmove((complexf*)input_buffer,((complexf*)input_buffer)+state.input_processed,(the_bufsize-state.input_processed)*sizeof(complexf)); //memmove lets the source and destination overlap
            fread(((complexf*)input_buffer)+(the_bufsize-state.input_processed), sizeof(complexf), state.input_processed, stdin);
            //fprintf(stderr,"iskip=%d state.output_size=%d start=%x target=%x skipcount=%x \n",state.input_processed,state.output_size,input_buffer, ((complexf*)input_buffer)+(BIG_BUFSIZE-state.input_processed),(BIG_BUFSIZE-state.input_processed));
        }
    }

    if(!strcmp(argv[1],"octave_complex_c"))
    {
        if(argc<=2) return badsyntax("need required parameter (samples_to_plot)");
        int samples_to_plot = 0;
        sscanf(argv[2], "%d", &samples_to_plot);
        if(samples_to_plot<=0) return badsyntax("Number of samples to plot should be > 0");
        if(argc<=3) return badsyntax("need required parameter (out_of_n_samples)");
        int out_of_n_samples = 0;
        sscanf(argv[3], "%d", &out_of_n_samples);
        if(out_of_n_samples<samples_to_plot) return badsyntax("out_of_n_samples should be < samples_to_plot");
        int mode2d = 0;
        if(argc>4) mode2d = !strcmp(argv[4], "--2d"); 
        complexf* read_buf = (complexf*)malloc(sizeof(complexf)*the_bufsize);

        if(!sendbufsize(initialize_buffers())) return -2;
        for(;;)
        {
            fread(read_buf, sizeof(complexf), samples_to_plot, stdin);          
            printf("N = %d;\nisig = [", samples_to_plot);
            for(int i=0;i<samples_to_plot;i++) printf("%f ", iof(read_buf, i));
            printf("];\nqsig = [");
            for(int i=0;i<samples_to_plot;i++) printf("%f ", qof(read_buf, i));
            printf("];\nzsig = [0:N-1];\n");
            if(mode2d) printf("subplot(2,1,1);\nplot(zsig,isig);\nsubplot(2,1,2);\nplot(zsig,qsig);\n");
            else printf("plot3(isig,zsig,qsig);\n");
            //printf("xlim([-1 1]);\nzlim([-1 1]);\n");
            fflush(stdout);
            //if(fseek(stdin, (out_of_n_samples - samples_to_plot)*sizeof(complexf), SEEK_CUR)<0) { perror("fseek error"); return -3; } //this cannot be used on stdin
            for(int seek_remain=out_of_n_samples-samples_to_plot;seek_remain>0;seek_remain-=samples_to_plot)
            {
                fread(read_buf, sizeof(complexf), MIN_M(samples_to_plot,seek_remain), stdin);
            }
            FEOF_CHECK;
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"psk_modulator_u8_c")) //<n_psk>
    {
        int n_psk;
        if(argc<=2) return badsyntax("need required parameter (n_psk)");
        sscanf(argv[2],"%d",&n_psk);
        if(n_psk<=0 || n_psk>256) return badsyntax("n_psk should be between 1 and 256");

        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize);

        for(;;)
        {
            FEOF_CHECK;
            fread((unsigned char*)input_buffer, sizeof(unsigned char), the_bufsize, stdin);
            psk_modulator_u8_c((unsigned char*)input_buffer, (complexf*)output_buffer, the_bufsize, n_psk);
            FWRITE_C;
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"duplicate_samples_ntimes_u8_u8")) //<sample_size_bytes> <ntimes>
    {
        int sample_size_bytes = 0, ntimes = 0;
        if(argc<=2) return badsyntax("need required parameter (sample_size_bytes)");
        sscanf(argv[2],"%d",&sample_size_bytes);    
        if(sample_size_bytes<=0) return badsyntax("sample_size_bytes should be >0");
        if(argc<=3) return badsyntax("need required parameter (ntimes)");
        sscanf(argv[3],"%d",&ntimes);   
        if(ntimes<=0) return badsyntax("ntimes should be >0");
        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize*ntimes);
        unsigned char* local_input_buffer = (unsigned char*)malloc(sizeof(unsigned char)*the_bufsize*sample_size_bytes);
        unsigned char* local_output_buffer = (unsigned char*)malloc(sizeof(unsigned char)*the_bufsize*sample_size_bytes*ntimes);
        for(;;)
        {
            FEOF_CHECK;
            fread((void*)local_input_buffer, sizeof(unsigned char), the_bufsize*sample_size_bytes, stdin);
            duplicate_samples_ntimes_u8_u8(local_input_buffer, local_output_buffer, the_bufsize*sample_size_bytes, sample_size_bytes, ntimes);
            fwrite((void*)local_output_buffer, sizeof(unsigned char), the_bufsize*sample_size_bytes*ntimes, stdout);
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"psk31_interpolate_sine_cc")) //<interpolation>
    {
        int interpolation;
        if(argc<=2) return badsyntax("need required parameter (interpolation)");
        sscanf(argv[2],"%d",&interpolation);    
        if(interpolation<=0) return badsyntax("interpolation should be >0"); 
        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize*interpolation);
        complexf* local_output_buffer = (complexf*)malloc(sizeof(complexf)*the_bufsize*interpolation);
        complexf last_input;
        iof(&last_input,0) = 0;
        qof(&last_input,0) = 0;
        for(;;)
        {
            FEOF_CHECK;
            FREAD_C;
            last_input = psk31_interpolate_sine_cc((complexf*)input_buffer, local_output_buffer, the_bufsize, interpolation, last_input);
            fwrite((void*)local_output_buffer, sizeof(complexf), the_bufsize*interpolation, stdout);
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"pack_bits_1to8_u8_u8")) 
    {
        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize*8);
        unsigned char* local_input_buffer = (unsigned char*)malloc(sizeof(unsigned char)*the_bufsize);
        unsigned char* local_output_buffer = (unsigned char*)malloc(sizeof(unsigned char)*the_bufsize*8);
        for(;;)
        {
            FEOF_CHECK;
            fread((void*)local_input_buffer, sizeof(unsigned char), the_bufsize, stdin);
            pack_bits_1to8_u8_u8(local_input_buffer, local_output_buffer, the_bufsize);
            fwrite((void*)local_output_buffer, sizeof(unsigned char), the_bufsize*8, stdout);
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"pack_bits_8to1_u8_u8")) 
    {
        if(!initialize_buffers()) return -2;
        sendbufsize(1);
        char local_input_buffer[8];
        for(;;)
        {
            FEOF_CHECK;
            fread((void*)local_input_buffer, sizeof(unsigned char), 8, stdin);
            unsigned char c = pack_bits_8to1_u8_u8(local_input_buffer);
            fwrite(&c, sizeof(unsigned char), 1, stdout);
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"psk31_varicode_encoder_u8_u8")) 
    {
        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize*8);
        int output_max_size=the_bufsize*30;
        int output_size;
        int input_processed;
        unsigned char* local_input_buffer = (unsigned char*)malloc(sizeof(unsigned char)*the_bufsize);
        unsigned char* local_output_buffer = (unsigned char*)malloc(sizeof(unsigned char)*output_max_size);
        fread((void*)local_input_buffer, sizeof(unsigned char), the_bufsize, stdin);
        for(;;)
        {
            psk31_varicode_encoder_u8_u8(local_input_buffer, local_output_buffer, the_bufsize, output_max_size, &input_processed, &output_size);
            //fprintf(stderr, "os = %d\n", output_size);
            fwrite((void*)local_output_buffer, sizeof(unsigned char), output_size, stdout);
            FEOF_CHECK;
            memmove(local_input_buffer, local_input_buffer+input_processed, the_bufsize-input_processed); 
            fread(input_buffer+the_bufsize-input_processed, sizeof(unsigned char), input_processed, stdin);
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"dump_u8")) 
    {
        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize*3);
        unsigned char* local_input_buffer = (unsigned char*)malloc(sizeof(unsigned char)*the_bufsize);
        for(;;)
        {
            FEOF_CHECK;
            fread((void*)local_input_buffer, sizeof(unsigned char), the_bufsize, stdin);
            for(int i=0;i<the_bufsize;i++) printf("%02x ", local_input_buffer[i]);
            TRY_YIELD;
        }
    }

    int differential_codec_encode = 0;
    if( (differential_codec_encode = !strcmp(argv[1],"differential_encoder_u8_u8")) || (!strcmp(argv[1],"differential_decoder_u8_u8")) )
    {
        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize);
        unsigned char* local_input_buffer = (unsigned char*)malloc(sizeof(unsigned char)*the_bufsize);
        unsigned char* local_output_buffer = (unsigned char*)malloc(sizeof(unsigned char)*the_bufsize);
        unsigned char state = 0;
        for(;;)
        {
            FEOF_CHECK;
            fread((void*)local_input_buffer, sizeof(unsigned char), the_bufsize, stdin);
            state = differential_codec(local_input_buffer, local_output_buffer, the_bufsize, differential_codec_encode, state);
            fwrite((void*)local_output_buffer, sizeof(unsigned char), the_bufsize, stdout);
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"bpsk_costas_loop_cc")) //<loop_bandwidth> <damping_factor> [--dd | --decision_directed] [--output_error | --output_dphase | --output_nco | --output_combined <error_file> <dphase_file> <nco_file>]
    {
        float loop_bandwidth;
        if(argc<=2) return badsyntax("need required parameter (loop_bandwidth)");
        sscanf(argv[2],"%f",&loop_bandwidth);

        float damping_factor;
        if(argc<=3) return badsyntax("need required parameter (damping_factor)");
        sscanf(argv[3],"%f",&damping_factor);

        int decision_directed = !!(argc>4 && (!strcmp(argv[4], "--dd") || !strcmp(argv[5], "--decision_directed")));
        if(decision_directed) { errhead(); fprintf(stderr, "decision directed mode\n"); }

        int output_error    =                  !!(argc>4+decision_directed && (!strcmp(argv[4+decision_directed], "--output_error")));
        int output_dphase   = !output_error  & !!(argc>4+decision_directed && (!strcmp(argv[4+decision_directed], "--output_dphase")));
        int output_nco      = !output_dphase & !!(argc>4+decision_directed && (!strcmp(argv[4+decision_directed], "--output_nco")));
        int output_combined = !output_nco    & !!(argc>4+decision_directed && (!strcmp(argv[4+decision_directed], "--output_combined")));


        bpsk_costas_loop_state_t state;
        init_bpsk_costas_loop_cc(&state, decision_directed, damping_factor, loop_bandwidth);
        errhead(); fprintf(stderr, "alpha = %f, beta = %f\n", state.alpha, state.beta);

        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize);

        float* buffer_output_error  = (!(output_combined || output_error))  ? NULL : (float*)malloc(sizeof(float)*the_bufsize);
        float* buffer_output_dphase = (!(output_combined || output_dphase)) ? NULL : (float*)malloc(sizeof(float)*the_bufsize);
        complexf* buffer_output_nco = (!(output_combined || output_nco))    ? NULL : (complexf*)malloc(sizeof(complexf)*the_bufsize);

        FILE* file_output_error = NULL;
        FILE* file_output_dphase = NULL;
        FILE* file_output_nco = NULL;
        if(output_combined)
        {
            if(!(argc>4+decision_directed+3)) { return badsyntax("need required parameters after --output_combined: <error_file> <dphase_file> <nco_file>"); }
            file_output_error  = fopen(argv[4+decision_directed+1], "w");
            file_output_dphase = fopen(argv[4+decision_directed+2], "w");
            file_output_nco    = fopen(argv[4+decision_directed+3], "w");
        }

        for(;;)
        {
            FEOF_CHECK;
            FREAD_C;
            bpsk_costas_loop_cc((complexf*)input_buffer, (complexf*)output_buffer, the_bufsize, 
                    buffer_output_error, buffer_output_dphase, buffer_output_nco, 
                    &state);
            if(output_error) fwrite(buffer_output_error, sizeof(float), the_bufsize, stdout);
            else if(output_dphase) fwrite(buffer_output_dphase, sizeof(float), the_bufsize, stdout);
            else if(output_nco) fwrite(buffer_output_nco, sizeof(complexf), the_bufsize, stdout);
            else 
            {
                if(output_combined) 
                {
                    fwrite(buffer_output_error, sizeof(float), the_bufsize, file_output_error);
                    fwrite(buffer_output_dphase, sizeof(float), the_bufsize, file_output_dphase);
                    fwrite(buffer_output_nco, sizeof(complexf), the_bufsize, file_output_nco);
                }
                FWRITE_C; 
            }
            TRY_YIELD;
        }
        fclose(file_output_error);
        fclose(file_output_dphase);
        fclose(file_output_nco);
    }

    if(!strcmp(argv[1],"simple_agc_cc")) //<rate> [reference [max_gain]] 
    {
        float rate;
        if(argc<=2) return badsyntax("need required parameter (rate)");
        sscanf(argv[2],"%f",&rate);
        if(rate<=0) return badsyntax("rate should be > 0");

        float reference = 1.;
        if(argc>3) sscanf(argv[3],"%f",&reference);
        if(reference<=0) return badsyntax("reference should be > 0");

        float max_gain = 65535.;
        if(argc>4) sscanf(argv[4],"%f",&max_gain);
        if(max_gain<=0) return badsyntax("max_gain should be > 0");

        float current_gain = 1.;

        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize);

        for(;;)
        {
            FEOF_CHECK;
            FREAD_C;
            simple_agc_cc((complexf*)input_buffer, (complexf*)output_buffer, the_bufsize, rate, reference, max_gain, &current_gain);
            FWRITE_C;
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"firdes_peak_c")) //<rate> <length> [window [--octave]]
    {
        //Process the params
        if(argc<=3) return badsyntax("need required parameters (rate, length)");

        float rate;
        sscanf(argv[2],"%g",&rate);
        int length;
        sscanf(argv[3],"%d",&length);
        if(length%2==0) return badsyntax("number of symmetric FIR filter taps should be odd");

        window_t window = WINDOW_DEFAULT;
        if(argc>=5)
        {
            window=firdes_get_window_from_string(argv[4]);
        }
        else { errhead(); fprintf(stderr,"window = %s\n",firdes_get_string_from_window(window)); }

        int octave=(argc>=6 && !strcmp("--octave",argv[5]));

        complexf* taps=(complexf*)malloc(sizeof(complexf)*length);

        //Make the filter
        firdes_add_peak_c(taps, length, rate, window, 0, 1);

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
            "subplot(2,1,2);plot(arg(fser));\n",fft_length-length);

        //Wait forever, so that octave won't close just after popping up the window.
        //You can close it with ^C.
        if(octave) { fflush(stdout); getchar(); }
        return 0;
    }
 
    if(!strcmp(argv[1],"peaks_fir_cc")) //<taps_length> <peak_rate × N>
    {
        //rule of thumb: bw = 2/taps_length,   which does not equal to transition_bw

        if(argc<=2) return badsyntax("need required parameter (taps_length)");
        int taps_length;
        sscanf(argv[2],"%d",&taps_length);

        int num_peaks = argc-3;
        float* peak_rate = (float*)malloc(sizeof(float)*num_peaks);
        for(int i=0;i<num_peaks;i++)
            sscanf(argv[3+i], "%f", peak_rate+i);
        if(num_peaks<=0) return badsyntax("need required parameter (peak_rate) once or multiple times");
        //for(int i=0;i<num_peaks;i++) fprintf(stderr, "%f\n", peak_rate[i]);
        fflush(stderr);

        window_t window = WINDOW_DEFAULT;

        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize);
        if(the_bufsize - taps_length <= 0 ) return badsyntax("taps_length is below buffer size, decrease taps_length");

        complexf* taps = (complexf*)calloc(sizeof(complexf),taps_length);
        for(int i=0; i<num_peaks; i++)
        {
            //fprintf(stderr, "nr = %d\n", i==num_peaks-1);
            firdes_add_peak_c(taps, taps_length, peak_rate[i], window, 1, i==num_peaks-1);
        }

        int output_size=0;
        FREAD_C;
        for(;;)
        {
            FEOF_CHECK;
            output_size = apply_fir_cc((complexf*)input_buffer, (complexf*)output_buffer, the_bufsize, taps, taps_length);
            fwrite(output_buffer, sizeof(complexf), output_size, stdout);
            //fprintf(stderr, "os = %d, is = %d\n", output_size, the_bufsize);
            TRY_YIELD;
            memmove((complexf*)input_buffer,((complexf*)input_buffer)+output_size,(the_bufsize-output_size)*sizeof(complexf)); 
            fread(((complexf*)input_buffer)+(the_bufsize-output_size), sizeof(complexf), output_size, stdin);
        }
    }

    if(!strcmp(argv[1], "repeat_u8"))
    {
        if(argc<=2) return badsyntax("no data to repeat");
        unsigned char* repeat_buffer = (unsigned char*)malloc(sizeof(unsigned char)*(argc-2));
        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize); //this is really (c-2) but this is a very fast source block so it makes no sense to send out a small number here
        for(int i=0;i<argc-2;i++)
        {
            FEOF_CHECK;
            int current_val;
            sscanf(argv[i+2], "%d", &current_val); 
            repeat_buffer[i]=current_val;
            TRY_YIELD;
        }
        for(;;) fwrite(repeat_buffer, sizeof(unsigned char), argc-2, stdout);
    }

    if(!strcmp(argv[1], "awgn_cc"))
    {
        FILE* urandom = init_get_random_samples_f();
        if(argc<=2) return badsyntax("required parameter <snr_db> is missing.");
        float snr_db = 0;
        sscanf(argv[2],"%f",&snr_db);
        FILE* awgnfile = NULL; 
        if(argc>=5 && !strcmp(argv[3],"--awgnfile")) 
        {
            awgnfile=fopen(argv[4], "r");
            if(!awgnfile) return badsyntax("failed to open the --awgnfile");
        }
        int parnumadd=2*(!!awgnfile);
        int snrshow = 0;
        if(argc>=4+parnumadd && !strcmp(argv[3+parnumadd],"--snrshow")) snrshow = 1;
        float signal_amplitude_per_noise = pow(10,snr_db/20);
        float a_signal=signal_amplitude_per_noise/(signal_amplitude_per_noise+1.0);
        float a_noise=1.0/(signal_amplitude_per_noise+1.0);
        errhead(); fprintf(stderr, "a_signal = %f, a_noise = %f\n", a_signal, a_noise);
        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize); 
        complexf* awgn_buffer = (complexf*)malloc(sizeof(complexf)*the_bufsize);
        for(;;)
        {
            FEOF_CHECK;
            FREAD_C;
            //get_awgn_samples_f((float*)awgn_buffer, the_bufsize*2, urandom);
            if(!awgnfile) get_random_gaussian_samples_c(awgn_buffer, the_bufsize, urandom);
            else 
            {
                for(;;)
                {
                    int items_read=fread(awgn_buffer, sizeof(complexf), the_bufsize, awgnfile);
                    if(items_read<the_bufsize) rewind(awgnfile);
                    else break;
                }
            }
            /*if(snrshow) 
            {
                float power_signal = total_logpower_cf((complexf*)input_buffer, the_bufsize);
                float power_noise = total_logpower_cf(awgn_buffer, the_bufsize);
                fprintf(stderr, "csdr awgn_cc: at the beginning, power_signal = %f dB, power_noise = %f dB\n", power_signal, power_noise);
            }*/
            gain_ff(input_buffer, input_buffer, the_bufsize*2, a_signal);
            gain_ff((float*)awgn_buffer, (float*)awgn_buffer, the_bufsize*2, a_noise*0.707);
            if(snrshow)
            {
                float power_signal = total_logpower_cf((complexf*)input_buffer, the_bufsize);
                float power_noise = total_logpower_cf(awgn_buffer, the_bufsize);
                //fprintf(stderr, "csdr awgn_cc: after gain_ff, power_signal = %f dB, power_noise = %f dB\n", power_signal, power_noise);
                errhead(); fprintf(stderr, "SNR = %f dB\n", power_signal - power_noise);
            }
            add_ff(input_buffer, (float*)awgn_buffer, (float*)output_buffer, the_bufsize*2);
            FWRITE_C;
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1], "uniform_noise_f"))
    {
        FILE* urandom = init_get_random_samples_f();
        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize); 
        for(;;)
        {
            FEOF_CHECK;
            get_random_samples_f(output_buffer, the_bufsize, urandom);
            FWRITE_R;
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1], "gaussian_noise_c"))
    {
        FILE* urandom = init_get_random_samples_f();
        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize); 
        for(;;)
        {
            FEOF_CHECK;
            get_random_gaussian_samples_c((complexf*)output_buffer, the_bufsize, urandom);
            FWRITE_C;
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1], "normalized_timing_variance_u32_f")) //<samples_per_symbol> <initial_sample_offset> [--debug]
    {
        int samples_per_symbol = 0;
        if(argc<=2) return badsyntax("required parameter <samples_per_symbol> is missing.");
        sscanf(argv[2],"%d",&samples_per_symbol);

        int initial_sample_offset = 0;
        if(argc<=3) return badsyntax("required parameter <initial_sample_offset> is missing.");
        sscanf(argv[3],"%d",&initial_sample_offset);

        int debug_print = 0;
        if(argc>4 && !strcmp(argv[4],"--debug")) debug_print = 1; 

        if(!initialize_buffers()) return -2;
        sendbufsize(the_bufsize); 
        float* temp_buffer = (float*)malloc(sizeof(float)*the_bufsize);
        for(;;)
        {
            FEOF_CHECK;
            FREAD_R; //doesn't count, reads 4 bytes per sample anyway
            float nv = normalized_timing_variance_u32_f((unsigned*)input_buffer, temp_buffer, the_bufsize, samples_per_symbol, initial_sample_offset, debug_print);
            fwrite(&nv, sizeof(float), 1, stdout);
            errhead(); fprintf(stderr, "normalized variance = %f\n", nv);
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1], "add_n_zero_samples_at_beginning_f")) //<n_zero_samples>
    {
        int n_zero_samples = 0;
        if(argc<=2) return badsyntax("required parameter <n_zero_samples> is missing.");
        sscanf(argv[2],"%d",&n_zero_samples);
        if(!sendbufsize(initialize_buffers())) return -2;
        float* zeros=(float*)calloc(sizeof(float),n_zero_samples);
        fwrite(zeros, sizeof(float), n_zero_samples, stdout);
        clone_(the_bufsize);
    }

    int pulse_shaping_filter_which = 0;
    if(
            (!strcmp(argv[1], "firdes_pulse_shaping_filter_f") && (pulse_shaping_filter_which = 1)) ||
            (!strcmp(argv[1], "pulse_shaping_filter_cc") && (pulse_shaping_filter_which = 2))
    ) //(RRC <samples_per_symbol> <num_taps> <beta> | COSINE <samples_per_symbol>)
    {
        if(argc<=2) return badsyntax("required parameter <pulse_shaping_filter_type> is missing.");
        matched_filter_type_t type = matched_filter_get_type_from_string(argv[2]);

        int samples_per_symbol = 0;
        if(argc<=3) return badsyntax("required parameter <samples_per_symbol> is missing.");
        sscanf(argv[3],"%d",&samples_per_symbol);

        int num_taps = 0;
        if(type!=MATCHED_FILTER_COSINE)
        {
            if(argc<=4) return badsyntax("required parameter <num_taps> is missing.");
            sscanf(argv[4],"%d",&num_taps);
        }
        else num_taps = (2*samples_per_symbol)+1;

        float beta = 0;
        if(type==MATCHED_FILTER_RRC)
        {
            if(argc<=5) return badsyntax("required parameter <beta> is missing.");
            sscanf(argv[5],"%f",&beta);
        }

        float* taps = (float*)malloc(sizeof(float)*num_taps);
        switch(type)
        {
        case MATCHED_FILTER_RRC:
            firdes_rrc_f(taps, num_taps, samples_per_symbol, beta);
            break;
        case MATCHED_FILTER_COSINE:
            firdes_cosine_f(taps, num_taps, samples_per_symbol);
            break;
        }
        //fprintf(stderr, "beta = %f, num_taps = %d, samples_per_symbol = %d\n", beta, num_taps, samples_per_symbol);

        if(!sendbufsize(initialize_buffers())) return -2;

        if(pulse_shaping_filter_which==1)
        {
            for(int i=0;i<num_taps;i++) printf("%f ", taps[i]);
            return 0;
        }

        int output_size=0;
        FREAD_C;
        for(;;)
        {
            FEOF_CHECK;
            output_size = apply_real_fir_cc((complexf*)input_buffer, (complexf*)output_buffer, the_bufsize, taps, num_taps);
            fwrite(output_buffer, sizeof(complexf), output_size, stdout);
            //fprintf(stderr, "os = %d, is = %d, num_taps = %d\n", output_size, the_bufsize, num_taps);
            TRY_YIELD;
            memmove((complexf*)input_buffer,((complexf*)input_buffer)+output_size,(the_bufsize-output_size)*sizeof(complexf)); 
            fread(((complexf*)input_buffer)+(the_bufsize-output_size), sizeof(complexf), output_size, stdin);
        }
    }

    if(!strcmp(argv[1], "generic_slicer_f_u8")) //<n_symbols>
    {
        int n_symbols = 0;
        if(argc<=2) return badsyntax("required parameter <n_symbols> is missing.");
        sscanf(argv[2],"%d",&n_symbols);
        if(!sendbufsize(initialize_buffers())) return -2;
        for(;;)
        {
            FEOF_CHECK;
            if(!FREAD_R) break;
            generic_slicer_f_u8(input_buffer, (unsigned char*)output_buffer, the_bufsize, n_symbols);
            FWRITE_U8;
            TRY_YIELD;
        }
        return 0;
    }

    if(!strcmp(argv[1], "plain_interpolate_cc")) //<interpolation>
    {
        int interpolation = 0;
        if(argc<=2) return badsyntax("required parameter <interpolation> is missing.");
        sscanf(argv[2],"%d",&interpolation);
        if(!sendbufsize(interpolation*initialize_buffers())) return -2;
        complexf* plainint_output_buffer = (complexf*)malloc(sizeof(complexf)*the_bufsize*interpolation);
        for(;;)
        {
            FEOF_CHECK;
            FREAD_C;
            plain_interpolate_cc((complexf*)input_buffer, plainint_output_buffer, the_bufsize, interpolation);
            fwrite(plainint_output_buffer, sizeof(float)*2, the_bufsize*interpolation, stdout);
            TRY_YIELD;
        }
        return 0;
    }

    if(!strcmp(argv[1], "dbpsk_decoder_c_u8")) 
    {
        if(!sendbufsize(initialize_buffers())) return -2;
        unsigned char* local_output_buffer = (unsigned char*)malloc(sizeof(unsigned char)*the_bufsize);
        for(;;)
        {
            FEOF_CHECK;
            FREAD_C;
            dbpsk_decoder_c_u8((complexf*)input_buffer, local_output_buffer, the_bufsize);
            fwrite(local_output_buffer, sizeof(unsigned char), the_bufsize, stdout);
            TRY_YIELD;
        }
        return 0;
    }

    if(!strcmp(argv[1], "bfsk_demod_cf")) //<spacing> <filter_length> 
    {
        float frequency_shift = 0;
        if(argc<=2) return badsyntax("required parameter <frequency_shift> is missing.");
        sscanf(argv[2],"%f",&frequency_shift);

        int filter_length = 0;
        if(argc<=3) return badsyntax("required parameter <filter_length> is missing.");
        sscanf(argv[3],"%d",&filter_length);

        complexf* mark_filter = (complexf*)malloc(sizeof(complexf)*filter_length);
        complexf* space_filter = (complexf*)malloc(sizeof(complexf)*filter_length);
        firdes_add_peak_c(mark_filter, filter_length, frequency_shift/2, WINDOW_DEFAULT, 0, 1);
        firdes_add_peak_c(space_filter, filter_length, -frequency_shift/2, WINDOW_DEFAULT, 0, 1);

        if(!sendbufsize(initialize_buffers())) return -2;

        int input_skip=0;
        int output_size=0;
        FREAD_C;
        for(;;)
        {
            FEOF_CHECK;
            output_size=bfsk_demod_cf((complexf*)input_buffer, output_buffer, the_bufsize, mark_filter, space_filter, filter_length);
            fwrite(output_buffer, sizeof(float), output_size, stdout);
            TRY_YIELD;
            memmove((complexf*)input_buffer,((complexf*)input_buffer)+output_size,(the_bufsize-output_size)*sizeof(complexf));
            fread(((complexf*)input_buffer)+(the_bufsize-output_size), sizeof(complexf), output_size, stdin);
        }
        return 0;
    }

    if(!strcmp(argv[1], "add_const_cc")) //<i> <q>
    {
        complexf x;
        if(argc<=2) return badsyntax("required parameter <add_i> is missing.");
        sscanf(argv[2],"%f",&iofv(x));
        if(argc<=2) return badsyntax("required parameter <add_q> is missing.");
        sscanf(argv[2],"%f",&qofv(x));

        if(!sendbufsize(initialize_buffers())) return -2;
        for(;;)
        {
            FEOF_CHECK;
            FREAD_C;
            add_const_cc((complexf*)input_buffer, (complexf*)output_buffer, the_bufsize, x);
            FWRITE_C;
            TRY_YIELD;
        }
        return 0;
    }

    if(!strcmp(argv[1], "tee")) //<path> [buffers]
    {
        if(argc<=2) return badsyntax("required parameter <path> is missing.");
        FILE* teefile = fopen(argv[2],"w");
        if(!teefile) return badsyntax("<path> cannot be opened!");
        errhead(); fprintf(stderr, "file opened: %s\n", argv[2]);
        int num_buffers=100;
        if(argc>3) sscanf(argv[3], "%d", &num_buffers);
        if(num_buffers<=0) return badsyntax("num_buffers should be > 0");
        SET_NONBLOCK(fileno(teefile));
        if(!sendbufsize(initialize_buffers())) return -2;
        unsigned char* async_tee_buffers = malloc(sizeof(unsigned char)*the_bufsize*num_buffers);
        int current_buffer_read_cntr = 0;
        int current_buffer_write_cntr = 0;
        int current_byte_write_cntr = 0;
        for(;;)
        {
            FEOF_CHECK;
            fread(async_tee_buffers+(the_bufsize*current_buffer_read_cntr), sizeof(unsigned char), the_bufsize, stdin);
            fwrite(async_tee_buffers+(the_bufsize*current_buffer_read_cntr++), sizeof(unsigned char), the_bufsize, stdout);
            if(current_buffer_read_cntr>=num_buffers) current_buffer_read_cntr = 0;
            if(current_buffer_read_cntr==current_buffer_write_cntr) { errhead(); fprintf(stderr, "circular buffer overflow (read pointer gone past write pointer)\n"); }
            //errhead(); fprintf(stderr, "new fwrites\n");
            while(current_buffer_write_cntr!=current_buffer_read_cntr)
            {
                int result = fwrite(async_tee_buffers+(the_bufsize*current_buffer_write_cntr)+current_byte_write_cntr, sizeof(unsigned char), the_bufsize-current_byte_write_cntr, teefile);
                if(!result) { errhead(); fprintf(stderr, "\t fwrite tee zero, next turn\n"); break; }
                current_byte_write_cntr += result;
                //errhead(); fprintf(stderr, "\tfwrite tee, current_byte_write_cntr = %d, current_buffer_write_cntr = %d, current_buffer_read_cntr = %d\n", 
                //        current_byte_write_cntr, current_buffer_write_cntr, current_buffer_read_cntr);
                if(current_byte_write_cntr >= the_bufsize) 
                {
                    current_byte_write_cntr = 0;
                    current_buffer_write_cntr++;
                    if(current_buffer_write_cntr>=num_buffers) current_buffer_write_cntr = 0;
                }
            }
            TRY_YIELD;
        }
        return 0;
    }

    if(!strcmp(argv[1],"shift_addition_fc"))
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
            errhead(); fprintf(stderr,"reinitialized to %g\n",rate); 
            int remain, current_size;
            float* ibufptr;
            float* obufptr;
            for(;;)
            {
                FEOF_CHECK;
                if(!FREAD_R) break;
                remain=the_bufsize;
                ibufptr=input_buffer;
                obufptr=output_buffer;
                while(remain)
                {
                    current_size=(remain>1024)?1024:remain;
                    starting_phase=shift_addition_fc(ibufptr, (complexf*)obufptr, current_size, data, starting_phase);
                    ibufptr+=current_size;
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

    if(!strcmp(argv[1],"fft_fc"))
    {
        /*
        For real FFT, the parameter is the number of output complex bins
        instead of the actual FFT size.
        Number of input samples used for each FFT is twice the given parameter.
        This makes it easier to replace fft_cc by fft_fc in some applications. */
        if(argc<=3) return badsyntax("need required parameters (fft_out_size, out_of_every_n_samples)");
        int fft_in_size=0, fft_out_size=0;
        sscanf(argv[2],"%d",&fft_out_size);
        if(log2n(fft_out_size)==-1) return badsyntax("fft_out_size should be power of 2");
        fft_in_size = 2*fft_out_size;
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
        sendbufsize(fft_out_size);

        //make FFT plan
        float* input=(float*)fft_malloc(sizeof(float)*fft_in_size);
        float* windowed=(float*)fft_malloc(sizeof(float)*fft_in_size);
        complexf* output=(complexf*)fft_malloc(sizeof(complexf)*fft_out_size);
        if(benchmark) { errhead(); fprintf(stderr,"benchmarking..."); }
        FFT_PLAN_T* plan=make_fft_r2c(fft_in_size, windowed, output, benchmark);
        if(benchmark) fprintf(stderr," done\n");
        //if(octave) printf("setenv(\"GNUTERM\",\"X11 noraise\");y=zeros(1,%d);semilogy(y,\"ydatasource\",\"y\");\n",fft_size); // TODO
        float *windowt;
        windowt = precalculate_window(fft_in_size, window);
        for(;;)
        {
            FEOF_CHECK;
            if(every_n_samples>fft_in_size)
            {
                fread(input, sizeof(float), fft_in_size, stdin);
                //skipping samples before next FFT (but fseek doesn't work for pipes)
                for(int seek_remain=every_n_samples-fft_in_size;seek_remain>0;seek_remain-=the_bufsize)
                {
                    fread(temp_f, sizeof(complexf), MIN_M(the_bufsize,seek_remain), stdin);
                }
            }
            else
            {
                //overlapped FFT
                for(int i=0;i<fft_in_size-every_n_samples;i++) input[i]=input[i+every_n_samples];
                fread(input+fft_in_size-every_n_samples, sizeof(float), every_n_samples, stdin);
            }
            //apply_window_c(input,windowed,fft_size,window);
            apply_precalculated_window_f(input,windowed,fft_in_size,windowt);
            fft_execute(plan);
            if(octave)
            {
#if 0
            // TODO
                printf("fftdata=[");
                //we have to swap the two parts of the array to get a valid spectrum
                for(int i=fft_size/2;i<fft_size;i++) printf("(%g)+(%g)*i ",iof(output,i),qof(output,i));
                for(int i=0;i<fft_size/2;i++) printf("(%g)+(%g)*i ",iof(output,i),qof(output,i));
                printf(
                    "];\n"
                    "y=abs(fftdata);\n"
                    "refreshdata;\n"
                );
#endif
            }
            else fwrite(output, sizeof(complexf), fft_out_size, stdout);
            TRY_YIELD;
        }
    }
/*
    if(!strcmp(argv[1],"syncword_search"))
    {
        if(argc<3) return badsyntax("need required parameter (syncword)");
        unsigned long syncword=0UL;
        int syncword_length=strlen(argv[2]);
        for(int i=0; i<sycword_len;i++)
        {
            syncword<<=4;
            char c=argv[2][i];
            unsigned char cval = 0;
            if(c<='9'&&c>='0') cval = c-'0';
            if(c<='f'&&c>='a') cval = c-'a'+10;
            syncword|=cval;
        }
        errhead(); fprintf("syncword = 0x%0x, syncword_length=%d\n", syncword, syncword_length);
        if(argc<4) return badsyntax("need required parameter (bits_after)");
        int bits_after = 0;
        sscanf(argv[3], &bits_after);
        if(bits_after<0) return badsyntax("bits_after should be >0");
        unsigned char* syncword_bits = malloc(sizeof(unsigned char)*syncword_length*4);
        int k=0;
        for(int i=0;i<syncword_length;i++)
        {
            for(int j=7;j;j--)
            {
                syncword_bits[k++]=syncword[i]>>j
            }
        }
        malloc

    }
*/
    if(!strcmp(argv[1],"pattern_search_u8_u8")) //<values_after> <pattern_values × N>
    {
        if(argc<3) return badsyntax("need required parameter (values_after)");
        int values_after = 0;
        sscanf(argv[2], "%d", &values_after);
        if(argc<4) return badsyntax("need required parameter (pattern_values × N)");
        int pattern_values_length = argc-3;
        unsigned* pattern_values = (unsigned*)malloc(sizeof(unsigned)*pattern_values_length);
        for(int i=0;i<pattern_values_length;i++)
            sscanf(argv[3+i],"%u",pattern_values+i);

        errhead(); fprintf(stderr,"pattern values: ");
        for(int i=0;i<pattern_values_length;i++)
            fprintf(stderr, "%x ", *(pattern_values+i));
        fprintf(stderr,"\n");

        unsigned char* input_buffer = (unsigned char*)malloc(sizeof(unsigned char)*pattern_values_length); //circular buffer
        unsigned char* output_buffer = (unsigned char*)malloc(sizeof(unsigned char)*values_after);
        int input_index = 0;
        int valid_values = 0;
        for(;;)
        {
            FEOF_CHECK;
            unsigned char cchar = input_buffer[input_index++]=(unsigned char)fgetc(stdin);
            if(valid_values<pattern_values_length) { valid_values++; continue; }
            if(input_index>=pattern_values_length) input_index=0;
            int match=1;
            //fprintf(stderr, "ov1: ");
            //for(int i=0;i<pattern_values_length;i++) fprintf(stderr, "%02x ",  input_buffer[i]);
            //fprintf(stderr, "\nov2: ");
            //for(int i=0;i<pattern_values_length;i++) fprintf(stderr, "%02x ",  pattern_values[i]);
            //fprintf(stderr, "\n");

            //fprintf(stderr, "v1: ");
            //for(int i=input_index;i<pattern_values_length;i++) fprintf(stderr, "%s%02x ", ((input_buffer[i])?"\x1B[34m":"\x1B[0m"), input_buffer[i]);
            //for(int i=0;i<input_index;i++) fprintf(stderr, "%s%02x ", ((input_buffer[i])?"\x1B[34m":"\x1B[0m"), input_buffer[i]);
            //fprintf(stderr, "\x1B[0m  %02x\n", cchar);

            //fprintf(stderr, "========\n");
            int j=0;
            for(int i=input_index;i<pattern_values_length;i++) 
            {
                //fprintf(stderr, "%02x ~ %02x\n",  input_buffer[i], pattern_values[j]);
                if(input_buffer[i]!=pattern_values[j++]) { match=0; break;}
            }
            //fprintf(stderr, "~~~~~~~~\n");
            if(input_index!=0 && match) 
            {
                for(int i=0;i<input_index;i++) 
                {
                    //fprintf(stderr, "%02x ~ %02x\n",  input_buffer[i], pattern_values[j]);
                    if(input_buffer[i]!=pattern_values[j++]) { match=0; break;}
                }
            }
           
            //if(match) fprintf(stderr, "j=%d\n", j);
            if(match ) 
            {
                valid_values = 0;
                //fprintf(stderr,"matched!\n");
                fread(output_buffer, sizeof(unsigned char), values_after, stdin);
                fwrite(output_buffer, sizeof(unsigned char), values_after, stdout);
            }
            TRY_YIELD;
        }
    }

    if(!strcmp(argv[1],"none"))
    {
        return 0;
    }

    if(argv[1][0]=='?' && argv[1][1]=='?')
    {
        char buffer[1000];
        snprintf(buffer, 1000-1, "xdg-open https://github.com/simonyiszk/csdr/blob/master/README.md#$(csdr ?%s | head -n1 | awk '{print $1;}')", argv[1]+2);
        fprintf(stderr, "csdr ??: %s\n", buffer);
        system(buffer);
        return 0;
    }

    if(argv[1][0]=='?')
    {
        char buffer[1000];
        snprintf(buffer, 1000-1, "csdr 2>&1 | grep -i %s", argv[1]+1);
        fprintf(stderr, "csdr ?: %s\n", buffer);
        system(buffer);
        return 0;
    }

    if(argv[1][0]=='=')
    {
        char buffer[100];
        snprintf(buffer, 100-1, "python -c \"import os, sys\nfrom math import *\nprint %s\"", argv[1]+1);
        system(buffer);
        return 0;
    }

    fprintf(stderr,"csdr: function name given in argument 1 (%s) does not exist. Possible causes:\n- You mistyped the commandline.\n- You need to update csdr to a newer version (if available).\n", argv[1]); return -1;
}
