/*
 This file is part of program wsprd, a detector/demodulator/decoder
 for the Weak Signal Propagation Reporter (WSPR) mode.
 
 Copyright 2001-2015, Joe Taylor, K1JT
 
 Much of the present code is based on work by Steven Franke, K9AN,
 which in turn was based on earlier work by K1JT.
 
 Copyright 2014-2015, Steven Franke, K9AN
 
 License: GNU GPL v3
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "types.h"
#include "misc.h"
#include "../libcsdr.h"
#include "fano.h"
#include "jelinek.h"

#include <assert.h>
#include <string.h>
#include <math.h>
#include <sched.h>
#include <time.h>
#include <fftw3.h>

#ifndef TRY_YIELD
 #define YIELD_EVERY_N_TIMES 64
 #define TRY_YIELD wspr_try_yield("TRY_YIELD", 0)
 #define TRY_YIELD_DECODE wspr_try_yield("TRY_YIELD_DECODE", 1)
#endif

// FIXME
#define ext_send_msg(...)
#define ext_send_msg_encoded(...)
#define ext_send_msg_data(...) 0

//#define WSPR_DEBUG_MSG	true
#define WSPR_DEBUG_MSG	false

//#define WSPR_PRINTF
#ifdef WSPR_PRINTF
	#define wprintf(fmt, ...) \
		fprintf(stderr, fmt, ## __VA_ARGS__)

	#define wdprintf(fmt, ...) \
	    fprintf(stderr, "WSPR %3ds ", timer_sec() - passes_start); \
		fprintf(stderr, fmt, ## __VA_ARGS__)
#else
	#define wprintf(fmt, ...)
	#define wdprintf(fmt, ...)
#endif

#define WSPR_FLOAT
#ifdef WSPR_FLOAT
	typedef float WSPR_REAL_t;
	typedef float WSPR_CPX_t;
    typedef struct {
        float re;
        float im;
    } WSPR_COMPLEX_t;

	#define WSPR_FFTW_COMPLEX fftwf_complex
	#define WSPR_FFTW_MALLOC fftwf_malloc
	#define WSPR_FFTW_FREE fftwf_free
	#define WSPR_FFTW_PLAN fftwf_plan
	#define WSPR_FFTW_PLAN_DFT_1D fftwf_plan_dft_1d
	#define WSPR_FFTW_DESTROY_PLAN fftwf_destroy_plan
	#define WSPR_FFTW_EXECUTE fftwf_execute
#else
	typedef double WSPR_REAL_t;
	typedef double WSPR_CPX_t;
    typedef struct {
        double re;
        double im;
    } WSPR_COMPLEX_t;

	#define WSPR_FFTW_COMPLEX fftw_complex
	#define WSPR_FFTW_MALLOC fftw_malloc
	#define WSPR_FFTW_FREE fftw_free
	#define WSPR_FFTW_PLAN fftw_plan
	#define WSPR_FFTW_PLAN_DFT_1D fftw_plan_dft_1d
	#define WSPR_FFTW_DESTROY_PLAN fftw_destroy_plan
	#define WSPR_FFTW_EXECUTE fftw_execute
#endif

#define K_PI (3.14159265358979323846)

#define	SYMTIME		(FSPS / FSRATE)		// symbol time: 256 samps @ 375 srate, 683 ms, 1.46 Hz

#define	SRATE		375					// design sample rate
#define	FSRATE		375.0
#define	CTIME		120					// capture time secs
#define	TPOINTS 	(SRATE * CTIME)

#define	FMIN		-110				// frequency range to search
#define	FMAX		110
//#define	FMIN		-150				// frequency range to search
//#define	FMAX		150
#define	BW_MAX		300.0				// +/- 150 Hz

// samples per symbol (at SRATE)
#define	FSPS		256.0				// round(SYMTIME * SRATE)
#define	SPS			256					// (int) FSPS
#define	NFFT		(SPS*2)
#define	HSPS		(SPS/2)

// groups
#define	GROUPS		(TPOINTS/NFFT)
#define	FPG			4					// FFTs per group

#define	NSYM_162	162
#define	FNSYM_162	162.0
#define	HSYM_81		(NSYM_162/2)
#define	FHSYM_81	(FNSYM_162/2)

#define	NBITS		HSYM_81

#define	LEN_DECODE	((NBITS + 7) / 8)
#define	LEN_CALL	(12 + SPACE_FOR_NULL)		// max of AANLL, ppp/AANLLL, AANLL/ss, plus 2
#define	LEN_C_L_P	(22 + SPACE_FOR_NULL)		// 22 = 12 call, sp, 6 grid, sp, 2 pwr
#define	LEN_GRID	(6 + SPACE_FOR_NULL)

#define	WSPR_STACKSIZE	2000

#define NPK 256
#define MAX_NPK 12
#define MAX_NPK_OLD 8

typedef struct {
	bool ignore;
	float freq0, snr0, drift0, sync0;
	int shift0, bin0;
	int freq_idx, flags;
	char snr_call[LEN_CALL];
} pk_t;

#define	WSPR_F_BIN			0x0fff
#define	WSPR_F_DECODING		0x1000
#define	WSPR_F_DELETE		0x2000
#define WSPR_F_DECODED		0x4000
#define	WSPR_F_IMAGE		0x8000

// assigned constants
extern int nffts;
extern int nbins_411;
extern int hbins_205;

typedef struct {
	float freq;
	char call[LEN_CALL];
	int hour, min;
	float snr, dt_print, drift1;
	double freq_print;
	char c_l_p[LEN_C_L_P];
} decode_t;

typedef struct {
	bool init;
	int rx_chan;
	int ping_pong, fft_ping_pong, decode_ping_pong;
	int capture;
	int status, status_resume;
	bool send_error, abort_decode;
	
	// csdr
	bool need_decode;
	float *input_buffer;
	int the_bufsize;
	FILE *demod_pipe;
	
	// options
	int quickmode, medium_effort, more_candidates, stackdecoder, subtraction;
	#define WSPR_TYPE_2MIN 2
	#define WSPR_TYPE_15MIN 15
	int wspr_type;
	
	// sampler
	bool reset, tsync;
	int last_min, last_sec;
	int decim, didx, group;
	double fi;
	
	// FFT task
	WSPR_FFTW_COMPLEX *fftin, *fftout;
	WSPR_FFTW_PLAN fftplan;
	int FFTtask_group;
	
	// computed by sampler or FFT task, processed by decode task
	#define N_PING_PONG 2
	time_t utc[N_PING_PONG];
	WSPR_CPX_t i_data[N_PING_PONG][TPOINTS], q_data[N_PING_PONG][TPOINTS];
	float pwr_samp[N_PING_PONG][NFFT][FPG*GROUPS];
	float pwr_sampavg[N_PING_PONG][NFFT];

	// decode task
	float min_snr, snr_scaling_factor;
	struct snode *stack;
	float dialfreq_MHz, cf_offset;
	u1_t symbols[NSYM_162], decdata[LEN_DECODE], channel_symbols[NSYM_162];
	char callsign[LEN_CALL], call_loc_pow[LEN_C_L_P], grid[LEN_GRID];
	decode_t deco[NPK];
} wspr_t;

// configuration
extern int bfo;

void wspr_try_yield(const char *from, int check_for_data);
void wspr_decode_init();
void wspr_init(char *demod_pipe_name, int decimate, int the_bufsize);
void wspr_data(WSPR_REAL_t *isamps, int nsamps);
void wspr_file(char *filename, int the_bufsize);
void wspr_decode_old(wspr_t *w);
void wspr_decode(wspr_t *w);
void wspr_send_peaks(wspr_t *w, pk_t *pk, int npk);
void wspr_hash_init();

void sync_and_demodulate(
	WSPR_CPX_t *id, WSPR_CPX_t *qd, long np,
	unsigned char *symbols, float *f1, int ifmin, int ifmax, float fstep,
	int *shift1,
	int lagmin, int lagmax, int lagstep,
	float drift1, int symfac, float *sync, int mode);

void renormalize(wspr_t *w, float psavg[], float smspec[]);

void unpack50(u1_t *dat, u4_t *call_28b, u4_t *grid_pwr_22b, u4_t *grid_15b, u4_t *pwr_7b);
int unpackcall(u4_t call_28b, char *call);
int unpackgrid(u4_t grid_15b, char *grid);
int unpackpfx(int32_t nprefix, char *call);
void deinterleave(unsigned char *sym);
int unpk_(u1_t *decdata, char *call_loc_pow, char *callsign, char *grid, int *dBm);
void subtract_signal(float *id, float *qd, long np,
	float f0, int shift0, float drift0, unsigned char* channel_symbols);
void subtract_signal2(float *id, float *qd, long np,
	float f0, int shift0, float drift0, unsigned char* channel_symbols);

int snr_comp(const void *elem1, const void *elem2);
int freq_comp(const void *elem1, const void *elem2);

typedef struct {
	double lat, lon;
} latLon_t;

void set_reporter_grid(char *grid);
double grid_to_distance_km(char *grid);
int latLon_to_grid6(latLon_t *loc, char *grid);
