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

#include "wspr.h"

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>

// wspr_status
#define	NONE		0
#define	IDLE		1
#define	SYNC		2
#define	RUNNING		3
#define	DECODING	4

static wspr_t wspr;

// assigned constants
int nffts, nbins_411, hbins_205;

// computed constants
static float window[NFFT];

// loaded from admin configuration
int bfo;

void wspr_try_yield(const char *from, int check_for_data)
{
    wspr_t *w = &wspr;

    if (check_for_data) {
        int n;
        do {
            // stdin was set non-blocking by csdr calling code
            n = fread(w->input_buffer, sizeof(float), w->the_bufsize, stdin);
            if (n > 0) {
                wspr_data(w->input_buffer, n/2);    // n/2 for number of complex samples
            } else {
                if (w->status != DECODING) {
                    // Not decoding: Sleep a little while so we don't hog the CPU.
                    // Otherwise the decoder needs to run at full speed and use regular
                    // Linux scheduling/time-slicing to share the CPU.
                    // This could be improved by not using read polling (FIXME).
                    usleep(100);
                }
            }
        } while (n > 0);
    }
    
    fflush(stdout);
    sched_yield();
}

const char *status_str[] = { "none", "idle", "sync", "running", "decoding" };

static void wspr_status(wspr_t *w, int status, int resume)
{
	wprintf("WSPR wspr_status: current %d-%s, set to %d-%s\n",
	    w->status, status_str[w->status], status, status_str[status]);
	
	// sends status as a control character (0x01 - 0x04) in the HTML stream
    fprintf(w->demod_pipe, "%c", status); fflush(w->demod_pipe);
	ext_send_msg(w->rx_chan, WSPR_DEBUG_MSG, "EXT WSPR_STATUS=%d", status);

	w->status = status;
	if (resume != NONE) {
		wprintf("WSPR wspr_status: will resume to %d-%s\n", resume, status_str[resume]);
		w->status_resume = resume;
	}
}

// compute FFTs incrementally during data capture rather than all at once at the beginning of the decode
void WSPR_FFT()
{
	int i,j,k;
	
    wspr_t *w = &wspr;
    int grp = w->FFTtask_group;
    int first = grp*FPG, last = first+FPG;
    //wprintf("WSPR FFT pp=%d grp=%d (%d-%d)\n", w->fft_ping_pong, grp, first, last);

    // Do ffts over 2 symbols, stepped by half symbols
    WSPR_CPX_t *id = w->i_data[w->fft_ping_pong], *qd = w->q_data[w->fft_ping_pong];

    float maxiq = 1e-66, maxpwr = 1e-66;
    int maxi=0;
    float savg[NFFT];
    memset(savg, 0, sizeof(savg));
    
    // FIXME: A large noise burst can washout the w->pwr_sampavg which prevents a proper
    // peak list being created in the decoder. An individual signal can be decoded fine
    // in the presence of the burst. But if there is no peak in the list the decoding
    // process is never started! We've seen this problem fairly frequently.
    
    for (i=first; i<last; i++) {
        for (j=0; j<NFFT; j++) {
            k = i*HSPS+j;
            w->fftin[j][0] = id[k] * window[j];
            w->fftin[j][1] = qd[k] * window[j];
            //if (i==0) wprintf("IN %d %fi %fq\n", j, w->fftin[j][0], w->fftin[j][1]);
        }
        
        //u4_t start = timer_us();
        TRY_YIELD;
        WSPR_FFTW_EXECUTE(w->fftplan);
        TRY_YIELD;
        //if (i==0) wprintf("512 FFT %.1f us\n", (float)(timer_us()-start));
        
        // NFFT = SPS*2
        // unwrap fftout:
        //      |1---> SPS
        // 2--->|      SPS

        for (j=0; j<NFFT; j++) {
            k = j+SPS;
            if (k > (NFFT-1))
                k -= NFFT;
            float ii = w->fftout[k][0];
            float qq = w->fftout[k][1];
            float pwr = ii*ii + qq*qq;
            //if (i==0) wprintf("OUT %d %fi %fq\n", j, ii, qq);
            if (ii > maxiq) { maxiq = ii; }
            if (pwr > maxpwr) { maxpwr = pwr; maxi = k; }
            w->pwr_samp[w->fft_ping_pong][j][i] = pwr;
            w->pwr_sampavg[w->fft_ping_pong][j] += pwr;
            savg[j] += pwr;
        }
        TRY_YIELD;
    }

    // send spectrum data to client
    float smspec[nbins_411];
    renormalize(w, savg, smspec);
    
    for (j=0; j<nbins_411; j++) {
        smspec[j] = 10*log10(smspec[j]) - w->snr_scaling_factor;
    }

    #define DATW nbins_411
    #define OUTW 1024
    #define LEVEL_COMP -30
    float out[OUTW];
    
    for (i = 0; i < OUTW; i++) {
        out[i] = smspec[(int) round((float) i / OUTW * (DATW-1))] + LEVEL_COMP;
    }

    fwrite(&out, sizeof(float), OUTW, stdout);

    //wprintf("WSPR_FFTtask group %d:%d %d-%d(%d) %f %f(%d)\n", w->fft_ping_pong, grp, first, last, nffts, maxiq, maxpwr, maxi);
    
    time_t t; time(&t); struct tm tm; gmtime_r(&t, &tm);
    wprintf("  %03d-%03d(%d) pp%d %02d:%02d\r", first, last, nffts, w->fft_ping_pong, tm.tm_min, tm.tm_sec); fflush(stderr);
}

void wspr_send_peaks(wspr_t *w, pk_t *pk, int npk)
{
    char peaks_s[NPK*(6+1 + LEN_CALL) + 16];
    char *s = peaks_s;
    *s = '\0';
    int j, n;
    for (j=0; j < npk; j++) {
    	int bin_flags = pk[j].bin0 | pk[j].flags;
    	n = sprintf(s, "%d:%s:", bin_flags, pk[j].snr_call); s += n;
    }
	ext_send_msg_encoded(w->rx_chan, WSPR_DEBUG_MSG, "EXT", "WSPR_PEAKS", "%s", peaks_s);
}

void WSPR_Decoder()
{
    u4_t start=0;

    wspr_t *w = &wspr;
    wspr_status(w, SYNC, SYNC);

    while (1) {
        if (w->need_decode) {
            w->need_decode = false;
            w->abort_decode = false;
        
            wspr_status(w, DECODING, NONE);
            start = timer_ms();
            wspr_decode(w);
            wprintf("WSPR Decoder TOTAL %.1f sec\n", (float)(timer_ms()-start)/1000.0);
        
            if (w->abort_decode)
                wprintf("WSPR decoder aborted\n");
        
            wspr_status(w, w->status_resume, NONE);
        }
        
        // important that this is the decode version that checks for input data
        TRY_YIELD_DECODE;
    }
}

static int int_decimate;

void wspr_data(WSPR_REAL_t *isamps, int nsamps)
{
	wspr_t *w = &wspr;
	int i;
    WSPR_COMPLEX_t *samps = (WSPR_COMPLEX_t *) isamps;
	
	//wprintf("WD%d didx %d send_error %d reset %d\n", w->capture, w->didx, w->send_error, w->reset);
	if (w->send_error) {
		wprintf("WSPR STOP send_error %d\n", w->send_error);
		w->send_error = FALSE;
		w->capture = FALSE;
		w->reset = TRUE;
	}

	if (w->reset) {
		w->ping_pong = w->decim = w->didx = w->group = 0;
		w->fi = 0;
		w->tsync = FALSE;
		w->status_resume = IDLE;	// decoder finishes after we stop capturing
		w->reset = FALSE;
	}

	if (!w->capture) {
		return;
	}
	
    // because of buffering this routine can be called with gaps > 1 sec (i.e. data gets bunched up)
	time_t t; time(&t); struct tm tm; gmtime_r(&t, &tm);

	if (tm.tm_sec != w->last_sec) {
        //fprintf(stderr, "WSPR %02d:%02d sync=%d\n", tm.tm_min, tm.tm_sec, w->tsync); fflush(stderr);
		if (tm.tm_min&1 && (tm.tm_sec >= 40 && tm.tm_sec <= 43))    // leave enough time for upload
			w->abort_decode = true;
	}
	
    if (w->tsync == FALSE) {		// sync to even minute boundary
        #define START_DELAY_SECS 1
        #define START_WINDOW 2
        if (!(tm.tm_min&1) && tm.tm_sec >= START_DELAY_SECS && tm.tm_sec <= (START_DELAY_SECS + START_WINDOW)) {

            if (w->status == RUNNING && !w->need_decode) {
                w->decode_ping_pong = w->ping_pong;
                w->need_decode = true;
            }

            w->ping_pong ^= 1;
            wprintf("\nWSPR SYNC ping_pong %d, %s", w->ping_pong, ctime(&t));
            w->decim = w->didx = w->group = 0;
            w->fi = 0;
            if (w->status != DECODING)
                wspr_status(w, RUNNING, RUNNING);
            w->tsync = TRUE;
        }
    }
	
	if (tm.tm_sec != w->last_sec) {
		w->last_min = tm.tm_min;
		w->last_sec = tm.tm_sec;
	}
	
	if (w->didx == 0) {
    	memset(&w->pwr_sampavg[w->ping_pong][0], 0, sizeof(w->pwr_sampavg[0]));
	}
	
	if (w->group == 0) w->utc[w->ping_pong] = t;
	
	WSPR_CPX_t *idat = w->i_data[w->ping_pong], *qdat = w->q_data[w->ping_pong];
	
    for (i=0; i<nsamps; i++) {

        // decimate
        if (w->decim++ < (int_decimate-1))
            continue;
        w->decim = 0;
        
        if (w->didx >= TPOINTS)
            return;

        if (w->group == 3) w->tsync = FALSE;	// arm re-sync
        
        WSPR_CPX_t re = (WSPR_CPX_t) samps[i].re;
        WSPR_CPX_t im = (WSPR_CPX_t) samps[i].im;
        idat[w->didx] = re;
        qdat[w->didx] = im;

        if ((w->didx % NFFT) == (NFFT-1)) {
            w->fft_ping_pong = w->ping_pong;
            w->FFTtask_group = w->group-1;
            if (w->group) WSPR_FFT();	// skip first to pipeline
            w->group++;
        }
        w->didx++;
    }
}

void wspr_init(char *demod_pipe_name, int decimate, int the_bufsize)
{
	int i;
    wspr_t *w;
    
    assert(FSPS == round(SYMTIME * FSRATE));
    assert(SPS == (int) FSPS);
    assert(HSPS == (SPS/2));

    nffts = FPG * floor(GROUPS-1) -1;
    nbins_411 = ceilf(NFFT * BW_MAX / FSRATE) +1;
    hbins_205 = (nbins_411-1)/2;

    wspr_decode_init();

    w = &wspr;
    memset(w, 0, sizeof(wspr_t));
    
    w->medium_effort = 1;
    w->wspr_type = WSPR_TYPE_2MIN;
    
    w->fftin = (WSPR_FFTW_COMPLEX*) WSPR_FFTW_MALLOC(sizeof(WSPR_FFTW_COMPLEX)*NFFT);
    w->fftout = (WSPR_FFTW_COMPLEX*) WSPR_FFTW_MALLOC(sizeof(WSPR_FFTW_COMPLEX)*NFFT);
    w->fftplan = WSPR_FFTW_PLAN_DFT_1D(NFFT, w->fftin, w->fftout, FFTW_FORWARD, FFTW_ESTIMATE);

    w->status_resume = IDLE;
    w->tsync = FALSE;
    w->capture = 0;
    time_t t; time(&t); struct tm tm; gmtime_r(&t, &tm);
    w->last_min = tm.tm_min;
    w->last_sec = tm.tm_sec;
    w->abort_decode = false;
    w->send_error = false;
	
	for (i=0; i < NFFT; i++) {
		window[i] = sin(i * K_PI/(NFFT-1));
	}

    // sdr source sampling rate determines if we rationally decimate here or rely on csdr fractional decimation
    int_decimate = decimate;
    wprintf("WSPR %s int_decimate=%d sps=%d NFFT=%d nbins_411=%d\n",
        (int_decimate != 1)? "LOCAL INTEGER DECIMATION" : "CSDR FRACTIONAL DECIMATION",
        int_decimate, SPS, NFFT, nbins_411);
    
    w->input_buffer = (float *) malloc(the_bufsize * sizeof(float));
    w->the_bufsize = the_bufsize;

    w->demod_pipe = fopen(demod_pipe_name, "w");
    if (w->demod_pipe == NULL) {
        wprintf("WSPR fopen FAIL %s\n", demod_pipe_name);
        exit(-1);
    }
    
    //fprintf(w->demod_pipe, "WSPR decoder running..\n"); fflush(w->demod_pipe);
    
    bfo = 750;
    w->dialfreq_MHz = 7040.1/1e3 - bfo/1e6;
    w->cf_offset = 100;
    wprintf("WSPR dialfreq=%f bfo=%d cf_offset=%f\n", w->dialfreq_MHz, bfo, w->cf_offset);
    w->rx_chan = 0;
    w->reset = true;
    w->capture = true;
	
	WSPR_Decoder();
}
