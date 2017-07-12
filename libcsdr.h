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
#define MAX_M(x,y) (((x)<(y))?(y):(x))

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
#define absof(complexf_input_p,i) (sqrt((iof(complexf_input_p,i)*iof(complexf_input_p,i))+(qof(complexf_input_p,i)*qof(complexf_input_p,i))))
#define argof(complexf_input_p,i) (atan2(qof(complexf_input_p,i),iof(complexf_input_p,i)))
#define cmult(cfo, cfi1, cfi2) {iof(cfo,0)=iof(cfi1,0)*iof(cfi2,0)-qof(cfi1,0)*qof(cfi2,0);qof(cfo,0)=iof(cfi1,0)*qof(cfi2,0)+iof(cfi2,0)*qof(cfi1,0);}
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

#define TIME_TAKEN(start,end) ((end.tv_sec-start.tv_sec)+(end.tv_nsec-start.tv_nsec)/1e9)

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
int firdes_filter_len(float transition_bw);

//demodulators
complexf fmdemod_quadri_cf(complexf* input, float* output, int input_size, float *temp, complexf last_sample);
complexf fmdemod_quadri_novect_cf(complexf* input, float* output, int input_size, complexf last_sample);
float fmdemod_atan_cf(complexf* input, float *output, int input_size, float last_phase);
void amdemod_cf(complexf* input, float *output, int input_size);
void amdemod_estimator_cf(complexf* input, float *output, int input_size, float alpha, float beta);
void limit_ff(float* input, float* output, int input_size, float max_amplitude);

//filters, decimators, resamplers, shift, etc.
float fir_one_pass_ff(float* input, float* taps, int taps_length);
int fir_decimate_cc(complexf *input, complexf *output, int input_size, int decimation, float *taps, int taps_length);
int fir_interpolate_cc(complexf *input, complexf *output, int input_size, int interpolation, float *taps, int taps_length);
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

float *precalculate_window(int size, window_t window);
void apply_window_c(complexf* input, complexf* output, int size, window_t window);
void apply_precalculated_window_c(complexf* input, complexf* output, int size, float *windowt);
void apply_precalculated_window_f(float* input, float* output, int size, float *windowt);
void apply_window_f(float* input, float* output, int size, window_t window);
void logpower_cf(complexf* input, float* output, int size, float add_db);
void accumulate_power_cf(complexf* input, float* output, int size);
void log_ff(float* input, float* output, int size, float add_db);

typedef struct fractional_decimator_ff_s
{
    float where;
    int input_processed;
    int output_size;
    int num_poly_points; //number of samples that the Lagrange interpolator will use
    float* poly_precalc_denomiator; //while we don't precalculate coefficients here as in a Farrow structure, because it is a fractional interpolator, but we rather precaculate part of the interpolator expression
    //float* last_inputs_circbuf; //circular buffer to store the last (num_poly_points) number of input samples.
    //int last_inputs_startsat; //where the circular buffer starts now
    //int last_inputs_samplewhere; 
    float* coeffs_buf;
    float* filtered_buf;
    int xifirst; 
    int xilast; 
    float rate;
    float *taps;
    int taps_length;
} fractional_decimator_ff_t;
fractional_decimator_ff_t fractional_decimator_ff_init(float rate, int num_poly_points, float* taps, int taps_length);
void fractional_decimator_ff(float* input, float* output, int input_size, fractional_decimator_ff_t* d);

typedef struct old_fractional_decimator_ff_s
{
    float remain;
    int input_processed;
    int output_size;
} old_fractional_decimator_ff_t;
old_fractional_decimator_ff_t old_fractional_decimator_ff(float* input, float* output, int input_size, float rate, float *taps, int taps_length, old_fractional_decimator_ff_t d);

typedef struct shift_table_data_s
{
    float* table;
    int table_size;
} shift_table_data_t;
void shift_table_deinit(shift_table_data_t table_data);
shift_table_data_t shift_table_init(int table_size);
float shift_table_cc(complexf* input, complexf* output, int input_size, float rate, shift_table_data_t table_data, float starting_phase);

typedef struct shift_addfast_data_s
{
    float dsin[4];
    float dcos[4];
    float phase_increment;
} shift_addfast_data_t;
shift_addfast_data_t shift_addfast_init(float rate);
shift_addfast_data_t shift_addfast_init(float rate);
float shift_addfast_cc(complexf *input, complexf* output, int input_size, shift_addfast_data_t* d, float starting_phase);

typedef struct shift_unroll_data_s
{
    float* dsin;
    float* dcos;
    float phase_increment;
    int size;
} shift_unroll_data_t;
float shift_unroll_cc(complexf *input, complexf* output, int input_size, shift_unroll_data_t* d, float starting_phase);
shift_unroll_data_t shift_unroll_init(float rate, int size);

int log2n(int x);
int next_pow2(int x);
void apply_fir_fft_cc(FFT_PLAN_T* plan, FFT_PLAN_T* plan_inverse, complexf* taps_fft, complexf* last_overlap, int overlap_size);
void gain_ff(float* input, float* output, int input_size, float gain);
float get_power_f(float* input, int input_size, int decimation);
float get_power_c(complexf* input, int input_size, int decimation);

void add_dcoffset_cc(complexf* input, complexf* output, int input_size);
float fmmod_fc(float* input, complexf* output, int input_size, float last_phase);
void fixed_amplitude_cc(complexf* input, complexf* output, int input_size, float amp);

void convert_u8_f(unsigned char* input, float* output, int input_size);
void convert_f_u8(float* input, unsigned char* output, int input_size);
void convert_s8_f(signed char* input, float* output, int input_size);
void convert_f_s8(float* input, signed char* output, int input_size);
void convert_f_s16(float* input, short* output, int input_size);
void convert_s16_f(short* input, float* output, int input_size);
void convert_f_i16(float* input, short* output, int input_size);
void convert_i16_f(short* input, float* output, int input_size);
void convert_f_s24(float* input, unsigned char* output, int input_size, int bigendian);
void convert_s24_f(unsigned char* input, float* output, int input_size, int bigendian);


int is_nan(float f);

//digital demod

typedef struct rtty_baudot_item_s
{
    unsigned long long code;
    unsigned char ascii_letter;
    unsigned char ascii_figure;
} rtty_baudot_item_t;

typedef enum rtty_baudot_decoder_state_e
{
    RTTY_BAUDOT_WAITING_STOP_PULSE = 0,
    RTTY_BAUDOT_WAITING_START_PULSE,
    RTTY_BAUDOT_RECEIVING_DATA
} rtty_baudot_decoder_state_t;

typedef struct rtty_baudot_decoder_s
{
    unsigned char fig_mode;
    unsigned char character_received;
    unsigned short shr;
    unsigned char bit_cntr;
    rtty_baudot_decoder_state_t state;
} rtty_baudot_decoder_t;

#define RTTY_FIGURE_MODE_SELECT_CODE 0b11011
#define RTTY_LETTER_MODE_SELECT_CODE 0b11111

char rtty_baudot_decoder_lookup(unsigned char* fig_mode, unsigned char c);
char rtty_baudot_decoder_push(rtty_baudot_decoder_t* s, unsigned char symbol);

//PSK31

typedef struct psk31_varicode_item_s
{
    unsigned long long code;
    int bitcount;
    unsigned char ascii;
} psk31_varicode_item_t;

char psk31_varicode_decoder_push(unsigned long long* status_shr, unsigned char symbol);

//Serial

typedef struct serial_line_s
{
    float samples_per_bits;
    int databits; //including parity
    float stopbits;
    int output_size;
    int input_used;
    float bit_sampling_width_ratio;
} serial_line_t;

void serial_line_decoder_f_u8(serial_line_t* s, float* input, unsigned char* output, int input_size);
void binary_slicer_f_u8(float* input, unsigned char* output, int input_size);


typedef enum pll_type_e
{
    PLL_P_CONTROLLER=1,
    PLL_PI_CONTROLLER=2
} pll_type_t;

typedef struct pll_s
{
    pll_type_t pll_type;
    //common:
    float output_phase;
    float dphase;
    float frequency;
    float alpha;
    float beta;
    float iir_temp;
} pll_t;

void pll_cc_init_pi_controller(pll_t* p, float bandwidth, float ko, float kd, float damping_factor);
void pll_cc_init_p_controller(pll_t* p, float alpha);
void pll_cc(pll_t* p, complexf* input, float* output_dphase, complexf* output_nco, int input_size);

typedef enum timing_recovery_algorithm_e
{
    TIMING_RECOVERY_ALGORITHM_GARDNER, 
    TIMING_RECOVERY_ALGORITHM_EARLYLATE 
} timing_recovery_algorithm_t;

#define TIMING_RECOVERY_ALGORITHM_DEFAULT TIMING_RECOVERY_ALGORITHM_GARDNER

typedef struct timing_recovery_state_s
{
    timing_recovery_algorithm_t algorithm;
    int decimation_rate; // = input_rate / output_rate. We should get an input signal that is N times oversampled. 
    int output_size;
    int input_processed;
    int use_q; //use both I and Q for calculating the error
    int debug_phase;
    int debug_every_nth;
    char* debug_writefiles_path;
    int last_correction_offset;
    float earlylate_ratio;
    float loop_gain;
    float max_error;
} timing_recovery_state_t;

timing_recovery_state_t timing_recovery_init(timing_recovery_algorithm_t algorithm, int decimation_rate, int use_q, float loop_gain, float max_error, int debug_every_nth, char* debug_writefiles_path);
void timing_recovery_cc(complexf* input, complexf* output, int input_size, float* timing_error, int* sampled_indexes,  timing_recovery_state_t* state);
timing_recovery_algorithm_t timing_recovery_get_algorithm_from_string(char* input);
char* timing_recovery_get_string_from_algorithm(timing_recovery_algorithm_t algorithm);
void octave_plot_point_on_cplxsig(complexf* signal, int signal_size, float error, int index, int correction_offset, char* writefiles_path, int points_size, ...);
void psk_modulator_u8_c(unsigned char* input, complexf* output, int input_size, int n_psk);
void duplicate_samples_ntimes_u8_u8(unsigned char* input, unsigned char* output, int input_size_bytes, int sample_size_bytes, int ntimes);
complexf psk31_interpolate_sine_cc(complexf* input, complexf* output, int input_size, int interpolation, complexf last_input);
void psk31_varicode_encoder_u8_u8(unsigned char* input, unsigned char* output, int input_size, int output_max_size, int* input_processed, int* output_size);
unsigned char differential_codec(unsigned char* input, unsigned char* output, int input_size, int encode, unsigned char state);

#if 0
typedef struct bpsk_costas_loop_state_s
{
    float rc_filter_alpha;
    float vco_phase_addition_multiplier;
    float vco_phase;
    float last_lpfi_output;
    float last_lpfq_output;
    float last_vco_phase_addition;
} bpsk_costas_loop_state_t;

bpsk_costas_loop_state_t init_bpsk_costas_loop_cc(float samples_per_bits);
void bpsk_costas_loop_cc(complexf* input, complexf* output, int input_size, bpsk_costas_loop_state_t* state);
#endif 

typedef struct bpsk_costas_loop_state_s
{
    float alpha;
    float beta;
    int decision_directed;
    float current_freq;
    float dphase;
    float nco_phase;
    float dphase_max;
    int dphase_max_reset_to_zero;
} bpsk_costas_loop_state_t;

void plain_interpolate_cc(complexf* input, complexf* output, int input_size, int interpolation);
void bpsk_costas_loop_cc(complexf* input, complexf* output, int input_size, float* output_error, float* output_dphase, complexf* output_nco, bpsk_costas_loop_state_t* s);
void init_bpsk_costas_loop_cc(bpsk_costas_loop_state_t* s, int decision_directed, float damping_factor, float bandwidth);

void simple_agc_cc(complexf* input, complexf* output, int input_size, float rate, float reference, float max_gain, float* current_gain);
void firdes_add_peak_c(complexf* output, int length, float rate, window_t window, int add, int normalize);
int apply_fir_cc(complexf* input, complexf* output, int input_size, complexf* taps, int taps_length);


FILE* init_get_random_samples_f();
void get_random_samples_f(float* output, int output_size, FILE* status);
void get_random_gaussian_samples_c(complexf* output, int output_size, FILE* status);
int deinit_get_random_samples_f(FILE* status);
float* add_ff(float* input1, float* input2, float* output, int input_size);
float total_logpower_cf(complexf* input, int input_size);
float normalized_timing_variance_u32_f(unsigned* input, float* temp, int input_size, int samples_per_symbol, int initial_sample_offset, int debug_print);

typedef enum matched_filter_type_e
{
    MATCHED_FILTER_RRC, 
    MATCHED_FILTER_COSINE 
} matched_filter_type_t;

#define MATCHED_FILTER_DEFAULT MATCHED_FILTER_RRC

int firdes_cosine_f(float* taps, int taps_length, int samples_per_symbol);
int firdes_rrc_f(float* taps, int taps_length, int samples_per_symbol, float beta);
matched_filter_type_t matched_filter_get_type_from_string(char* input);
int apply_real_fir_cc(complexf* input, complexf* output, int input_size, float* taps, int taps_length);
void generic_slicer_f_u8(float* input, unsigned char* output, int input_size, int n_symbols);
void plain_interpolate_cc(complexf* input, complexf* output, int input_size, int interpolation);;
void normalize_fir_f(float* input, float* output, int length);
float* add_const_cc(complexf* input, complexf* output, int input_size, complexf x);
void pack_bits_1to8_u8_u8(unsigned char* input, unsigned char* output, int input_size);
unsigned char pack_bits_8to1_u8_u8(unsigned char* input);
void dbpsk_decoder_c_u8(complexf* input, unsigned char* output, int input_size);
int bfsk_demod_cf(complexf* input, float* output, int input_size, complexf* mark_filter, complexf* space_filter, int taps_length);
