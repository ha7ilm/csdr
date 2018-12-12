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

static const unsigned char pr3[NSYM_162]=
{1,1,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,
    0,1,0,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,1,
    0,0,0,0,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,
    1,0,1,0,0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,1,
    0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,1,0,
    0,0,0,0,1,0,0,1,0,0,1,1,1,0,1,1,0,0,1,1,
    0,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,1,1,
    0,0,0,0,0,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,
    0,0};

void sync_and_demodulate(
	WSPR_CPX_t *id, WSPR_CPX_t *qd, long np,
	unsigned char *symbols, float *f1, int ifmin, int ifmax, float fstep,
	int *shift1,
	int lagmin, int lagmax, int lagstep,
	float drift1, int symfac, float *sync, int mode)
{
    /***********************************************************************
     * mode = 0: no frequency or drift search. find best time lag.          *
     *        1: no time lag or drift search. find best frequency.          *
     *        2: no frequency or time lag search. calculate soft-decision   *
     *           symbols using passed frequency and shift.                  *
     ************************************************************************/
    
    static float fplast=-10000.0;
    static float dt=1.0/FSRATE, df=FSRATE/FSPS;
    static float pi=K_PI;
    float twopidt, df15=df*1.5, df05=df*0.5;

    int i, j, k, lag;
    float
		i0[NSYM_162], q0[NSYM_162],
		i1[NSYM_162], q1[NSYM_162],
		i2[NSYM_162], q2[NSYM_162],
		i3[NSYM_162], q3[NSYM_162];
    double p0, p1, p2, p3, cmet, totp, syncmax, fac;
    double
		c0[SPS], s0[SPS],
		c1[SPS], s1[SPS],
		c2[SPS], s2[SPS],
		c3[SPS], s3[SPS];
    double
		dphi0, cdphi0, sdphi0,
		dphi1, cdphi1, sdphi1,
		dphi2, cdphi2, sdphi2,
		dphi3, cdphi3, sdphi3;
    float f0=0.0, fp, ss, fbest=0.0, fsum=0.0, f2sum=0.0, fsymb[NSYM_162];
    int best_shift = 0, ifreq;
    
    syncmax=-1e30;

    if( mode == 0 ) { ifmin=0; ifmax=0; fstep=0.0; f0=*f1; }
    if( mode == 1 ) { lagmin=*shift1; lagmax=*shift1; f0=*f1; }
    if( mode == 2 ) { lagmin=*shift1; lagmax=*shift1; ifmin=0; ifmax=0; f0=*f1; }

    twopidt=2*pi*dt;
    for(ifreq=ifmin; ifreq<=ifmax; ifreq++) {
        f0=*f1+ifreq*fstep;
        for(lag=lagmin; lag<=lagmax; lag=lag+lagstep) {
            ss=0.0;
            totp=0.0;
            for (i=0; i<NSYM_162; i++) {
                fp = f0 + (drift1/2.0)*((float)i-FHSYM_81)/FHSYM_81;
                if( i==0 || (fp != fplast) ) {  // only calculate sin/cos if necessary
                    dphi0=twopidt*(fp-df15);
                    cdphi0=cos(dphi0);
                    sdphi0=sin(dphi0);
                    
                    dphi1=twopidt*(fp-df05);
                    cdphi1=cos(dphi1);
                    sdphi1=sin(dphi1);
                    
                    dphi2=twopidt*(fp+df05);
                    cdphi2=cos(dphi2);
                    sdphi2=sin(dphi2);
                    
                    dphi3=twopidt*(fp+df15);
                    cdphi3=cos(dphi3);
                    sdphi3=sin(dphi3);
                    
                    c0[0]=1; s0[0]=0;
                    c1[0]=1; s1[0]=0;
                    c2[0]=1; s2[0]=0;
                    c3[0]=1; s3[0]=0;
                    
                    for (j=1; j<SPS; j++) {
                        c0[j]=c0[j-1]*cdphi0 - s0[j-1]*sdphi0;
                        s0[j]=c0[j-1]*sdphi0 + s0[j-1]*cdphi0;
                        c1[j]=c1[j-1]*cdphi1 - s1[j-1]*sdphi1;
                        s1[j]=c1[j-1]*sdphi1 + s1[j-1]*cdphi1;
                        c2[j]=c2[j-1]*cdphi2 - s2[j-1]*sdphi2;
                        s2[j]=c2[j-1]*sdphi2 + s2[j-1]*cdphi2;
                        c3[j]=c3[j-1]*cdphi3 - s3[j-1]*sdphi3;
                        s3[j]=c3[j-1]*sdphi3 + s3[j-1]*cdphi3;
						
						if ((j%YIELD_EVERY_N_TIMES)==(YIELD_EVERY_N_TIMES-1)) TRY_YIELD_DECODE;
                    }
                    fplast = fp;
                }
                
                i0[i]=0.0; q0[i]=0.0;
                i1[i]=0.0; q1[i]=0.0;
                i2[i]=0.0; q2[i]=0.0;
                i3[i]=0.0; q3[i]=0.0;
                
                for (j=0; j<SPS; j++) {
                    k=lag+i*SPS+j;
                    if( (k>0) && (k<np) ) {
                        i0[i]=i0[i] + id[k]*c0[j] + qd[k]*s0[j];
                        q0[i]=q0[i] - id[k]*s0[j] + qd[k]*c0[j];
                        i1[i]=i1[i] + id[k]*c1[j] + qd[k]*s1[j];
                        q1[i]=q1[i] - id[k]*s1[j] + qd[k]*c1[j];
                        i2[i]=i2[i] + id[k]*c2[j] + qd[k]*s2[j];
                        q2[i]=q2[i] - id[k]*s2[j] + qd[k]*c2[j];
                        i3[i]=i3[i] + id[k]*c3[j] + qd[k]*s3[j];
                        q3[i]=q3[i] - id[k]*s3[j] + qd[k]*c3[j];
						
						if ((j%YIELD_EVERY_N_TIMES)==(YIELD_EVERY_N_TIMES-1)) TRY_YIELD_DECODE;
                    }
                }
                p0=i0[i]*i0[i] + q0[i]*q0[i];
                p1=i1[i]*i1[i] + q1[i]*q1[i];
                p2=i2[i]*i2[i] + q2[i]*q2[i];
                p3=i3[i]*i3[i] + q3[i]*q3[i];

                p0=sqrt(p0);
                p1=sqrt(p1);
                p2=sqrt(p2);
                p3=sqrt(p3);
                
                totp=totp+p0+p1+p2+p3;
                cmet=(p1+p3)-(p0+p2);
                ss = (pr3[i] == 1) ? ss+cmet : ss-cmet;
                if( mode == 2) {                 //Compute soft symbols
                    if(pr3[i]==1) {
                        fsymb[i]=p3-p1;
                    } else {
                        fsymb[i]=p2-p0;
                    }
                }
            }
            ss=ss/totp;
            if( ss > syncmax ) {          //Save best parameters
                syncmax=ss;
                best_shift=lag;
                fbest=f0;
            }
        } // lag loop
    } //freq loop
    
    if( mode <=1 ) {                       //Send best params back to caller
        *sync=syncmax;
        *shift1=best_shift;
        *f1=fbest;
        return;
    }
    
    if( mode == 2 ) {
        *sync=syncmax;
        for (i=0; i<NSYM_162; i++) {              //Normalize the soft symbols
            fsum=fsum+fsymb[i]/FNSYM_162;
            f2sum=f2sum+fsymb[i]*fsymb[i]/FNSYM_162;
        }
        fac=sqrt(f2sum-fsum*fsum);
        TRY_YIELD_DECODE;

        for (i=0; i<NSYM_162; i++) {
            fsymb[i]=symfac*fsymb[i]/fac;
            if( fsymb[i] > 127) fsymb[i]=127.0;
            if( fsymb[i] < -128 ) fsymb[i]=-128.0;
            symbols[i] = fsymb[i] + 128;
        }
        TRY_YIELD_DECODE;
        return;
    }
    TRY_YIELD_DECODE;
    return;
}

// uses TRY_YIELD instead of TRY_YIELD_DECODE since called from non-decoder FFT code
void renormalize(wspr_t *w, float psavg[], float smspec[])
{
	int i,j,k;

	// smooth with 7-point window and limit the spectrum to +/-150 Hz
    int window[7] = {1,1,1,1,1,1,1};
    
    for (i=0; i<nbins_411; i++) {
        smspec[i] = 0.0;
        for (j=-3; j<=3; j++) {
            k = SPS-hbins_205+i+j;
            smspec[i] += window[j+3]*psavg[k];
        }
    }
	TRY_YIELD;

	// Sort spectrum values, then pick off noise level as a percentile
    float tmpsort[nbins_411];
    for (j=0; j<nbins_411; j++)
        tmpsort[j] = smspec[j];
    qsort(tmpsort, nbins_411, sizeof(float), qsort_floatcomp);
	TRY_YIELD;

	// Noise level of spectrum is estimated as 123/411= 30'th percentile
    float noise_level = tmpsort[122];
    
	/* Renormalize spectrum so that (large) peaks represent an estimate of snr.
	 * We know from experience that threshold snr is near -7dB in wspr bandwidth,
	 * corresponding to -7-26.3=-33.3dB in 2500 Hz bandwidth.
	 * The corresponding threshold is -42.3 dB in 2500 Hz bandwidth for WSPR-15.
	 */
	w->min_snr = pow(10.0,-7.0/10.0); //this is min snr in wspr bw
	w->snr_scaling_factor = (w->wspr_type == WSPR_TYPE_2MIN)? 26.3 : 35.3;
	for (j=0; j<411; j++) {
		smspec[j]=smspec[j]/noise_level - 1.0;
		if( smspec[j] < w->min_snr) smspec[j]=0.1*w->min_snr;
	}
	TRY_YIELD;
}

/***************************************************************************
 symbol-by-symbol signal subtraction
 ****************************************************************************/
void subtract_signal(float *id, float *qd, long np,
                     float f0, int shift0, float drift0, unsigned char* channel_symbols)
{
    float dt=1.0/FSRATE, df=FSRATE/FSPS;
    int i, j, k;
    float pi=K_PI, twopidt, fp;
    
    float i0,q0;
    float c0[SPS],s0[SPS];
    float dphi, cdphi, sdphi;
    
    twopidt=2*pi*dt;
    
    for (i=0; i<NSYM_162; i++) {
        fp = f0 + ((float)drift0/2.0)*((float)i-FHSYM_81)/FHSYM_81;
        
        dphi=twopidt*(fp+((float)channel_symbols[i]-1.5)*df);
        cdphi=cos(dphi);
        sdphi=sin(dphi);
        
        c0[0]=1; s0[0]=0;
        
        for (j=1; j<SPS; j++) {
            c0[j]=c0[j-1]*cdphi - s0[j-1]*sdphi;
            s0[j]=c0[j-1]*sdphi + s0[j-1]*cdphi;
        }
        
        i0=0.0; q0=0.0;
        
        for (j=0; j<SPS; j++) {
            k=shift0+i*SPS+j;
            if( (k>0) & (k<np) ) {
                i0=i0 + id[k]*c0[j] + qd[k]*s0[j];
                q0=q0 - id[k]*s0[j] + qd[k]*c0[j];
            }
        }
        
        
        // subtract the signal here.
        
        i0=i0/FSPS; //will be wrong for partial symbols at the edges...
        q0=q0/FSPS;
        
        for (j=0; j<SPS; j++) {
            k=shift0+i*SPS+j;
            if( (k>0) & (k<np) ) {
                id[k]=id[k]- (i0*c0[j] - q0*s0[j]);
                qd[k]=qd[k]- (q0*c0[j] + i0*s0[j]);
            }
        }
    }
    return;
}

/******************************************************************************
 Fully coherent signal subtraction
 *******************************************************************************/
void subtract_signal2(float *id, float *qd, long np,
                      float f0, int shift0, float drift0, unsigned char* channel_symbols)
{
    float dt=1.0/FSRATE, df=FSRATE/FSPS;
    float pi=K_PI, twopidt, phi=0, dphi, cs;
    int i, j, k, ii, nsym=NSYM_162, nspersym=SPS,  nfilt=SPS; //nfilt must be even number.
    int nsig=nsym*nspersym;
    int nc2=TPOINTS;
    
    float *refi, *refq, *ci, *cq, *cfi, *cfq;

    refi = (float *) malloc(sizeof(float)*nc2);
    refq = (float *) malloc(sizeof(float)*nc2);
    ci = (float *) malloc(sizeof(float)*nc2);
    cq = (float *) malloc(sizeof(float)*nc2);
    cfi = (float *) malloc(sizeof(float)*nc2);
    cfq = (float *) malloc(sizeof(float)*nc2);
    
    memset(refi,0,sizeof(float)*nc2);
    memset(refq,0,sizeof(float)*nc2);
    memset(ci,0,sizeof(float)*nc2);
    memset(cq,0,sizeof(float)*nc2);
    memset(cfi,0,sizeof(float)*nc2);
    memset(cfq,0,sizeof(float)*nc2);
    
    twopidt=2.0*pi*dt;
    
    /******************************************************************************
     Measured signal:                    s(t)=a(t)*exp( j*theta(t) )
     Reference is:                       r(t) = exp( j*phi(t) )
     Complex amplitude is estimated as:  c(t)=LPF[s(t)*conjugate(r(t))]
     so c(t) has phase angle theta-phi
     Multiply r(t) by c(t) and subtract from s(t), i.e. s'(t)=s(t)-c(t)r(t)
     *******************************************************************************/
    
    // create reference wspr signal vector, centered on f0.
    //
    for (i=0; i<nsym; i++) {
        
        cs=(float)channel_symbols[i];
        
        dphi=twopidt*
        (
         f0 + (drift0/2.0)*((float)i-(float)nsym/2.0)/((float)nsym/2.0)
         + (cs-1.5)*df
         );
        
        for ( j=0; j<nspersym; j++ ) {
            ii=nspersym*i+j;
            refi[ii]=cos(phi); //cannot precompute sin/cos because dphi is changing
            refq[ii]=sin(phi);
            phi=phi+dphi;
        }
    }
    
    // s(t) * conjugate(r(t))
    // beginning of first symbol in reference signal is at i=0
    // beginning of first symbol in received data is at shift0.
    // filter transient lasts nfilt samples
    // leave nfilt zeros as a pad at the beginning of the unfiltered reference signal
    for (i=0; i<nsym*nspersym; i++) {
        k=shift0+i;
        if( (k>0) && (k<np) ) {
            ci[i+nfilt] = id[k]*refi[i] + qd[k]*refq[i];
            cq[i+nfilt] = qd[k]*refi[i] - id[k]*refq[i];
        }
    }
    
    //lowpass filter and remove startup transient
    float w[nfilt], norm=0, partialsum[nfilt];
    memset(partialsum,0,sizeof(float)*nfilt);
    for (i=0; i<nfilt; i++) {
        w[i]=sin(pi*(float)i/(float)(nfilt-1));
        norm=norm+w[i];
    }
    for (i=0; i<nfilt; i++) {
        w[i]=w[i]/norm;
    }
    for (i=1; i<nfilt; i++) {
        partialsum[i]=partialsum[i-1]+w[i];
    }
    
    // LPF
    for (i=nfilt/2; i<TPOINTS-nfilt/2; i++) {
        cfi[i]=0.0; cfq[i]=0.0;
        for (j=0; j<nfilt; j++) {
            cfi[i]=cfi[i]+w[j]*ci[i-nfilt/2+j];
            cfq[i]=cfq[i]+w[j]*cq[i-nfilt/2+j];
        }
    }
    
    // subtract c(t)*r(t) here
    // (ci+j*cq)(refi+j*refq)=(ci*refi-cq*refq)+j(ci*refq)+cq*refi)
    // beginning of first symbol in reference signal is at i=nfilt
    // beginning of first symbol in received data is at shift0.
    for (i=0; i<nsig; i++) {
        if( i<nfilt/2 ) {        // take care of the end effect (LPF step response) here
            norm=partialsum[nfilt/2+i];
        } else if( i>(nsig-1-nfilt/2) ) {
            norm=partialsum[nfilt/2+nsig-1-i];
        } else {
            norm=1.0;
        }
        k=shift0+i;
        j=i+nfilt;
        if( (k>0) && (k<np) ) {
            id[k]=id[k] - (cfi[j]*refi[i]-cfq[j]*refq[i])/norm;
            qd[k]=qd[k] - (cfi[j]*refq[i]+cfq[j]*refi[i])/norm;
        }
    }
    
    free(refi);
    free(refq);
    free(ci);
    free(cq);
    free(cfi);
    free(cfq);

    return;
}

#include "./metric_tables.h"
static int mettab[2][256];
    
void wspr_decode_init()
{
    float bias=0.45;						//Fano metric bias (used for both Fano and stack algorithms)

    // setup metric table
    for (int i=0; i < SPS; i++) {
        mettab[0][i] = round(10 * (metric_tables[2][i] - bias));
        mettab[1][i] = round(10 * (metric_tables[2][SPS-1-i] - bias));
    }
    
    wspr_hash_init();
}
    
void wspr_decode(wspr_t *w)
{
    char cr[] = "(C) 2016, Steven Franke - K9AN";
    (void) cr;

    int i,j,k;

    int ipass, npasses = 1;
    int shift1, lagmin, lagmax, lagstep, ifmin, ifmax;
    unsigned int npoints = TPOINTS, metric, cycles, maxnp;

    float df = FSRATE/FSPS/2;
    float dt = 1.0/FSRATE, dt_print;

	int pki, npk;
    pk_t pk[NPK], pk_freq[NPK];

    double freq_print;
    float f1, fstep, sync1, drift1, snr;

    int ndecodes_pass;
	u4_t passes_start = timer_sec();
    
    //jksd FIXME need way to select:
    //	more_candidates
    //	stackdecoder
    //	subtraction
    
    // Parameters used for performance-tuning:
    unsigned int maxcycles=200;				//Decoder timeout limit
    float minsync1=0.10;					//First sync limit
    float minsync2=0.12;					//Second sync limit
    int iifac=2;							//Step size in final DT peakup
    int jig_range=128;
    int symfac=50;							//Soft-symbol normalizing factor
    int maxdrift=4;							//Maximum (+/-) drift
    float minrms=52.0 * (symfac/64.0);		//Final test for plausible decoding
    int delta=60;							//Fano threshold step
    
	wprintf("WSPR DECODE using decode_ping_pong %d\n", w->decode_ping_pong);

	WSPR_CPX_t *idat = w->i_data[w->decode_ping_pong];
	WSPR_CPX_t *qdat = w->q_data[w->decode_ping_pong];

    static bool wspr_decode_init;
    if (!wspr_decode_init) {
		if (w->stackdecoder)
			w->stack = (struct snode *) malloc(WSPR_STACKSIZE * sizeof(struct snode));

    	wspr_decode_init = true;
    }

	int uniques = 0;
	
	// multi-pass strategies
	//#define SUBTRACT_SIGNAL		// FIXME: how to implement spectrum subtraction given our incrementally-computed FFTs?
	#define MORE_EFFORT			// this scheme repeats work as maxcycles is increased, but it's difficult to eliminate that
		
	#if defined(SUBTRACT_SIGNAL)
		npasses = 2;
	#elif defined(MORE_EFFORT)
		npasses = 0;	// unlimited
	#endif

    for (ipass=0; (npasses == 0 || ipass < npasses) && !w->abort_decode; ipass++) {

		#if defined(SUBTRACT_SIGNAL)
        	if (ipass > 0 && ndecodes_pass == 0) break;
        #elif defined(MORE_EFFORT)
        	if (ipass == 0) {
        		maxcycles = 200;
        		iifac = 2;
        	} else {
        		if (maxcycles < 10000) maxcycles *= 2;
        		iifac = 1;
        	}
        #endif
        
		// only build the peak list on the first pass
		if (ipass == 0) {
			float smspec[nbins_411];
			renormalize(w, w->pwr_sampavg[w->decode_ping_pong], smspec);
			
			// Find all local maxima in smoothed spectrum.
			memset(pk, 0, sizeof(pk));
			
			npk = 0;
			bool candidate;
			
			if (w->more_candidates) {
				for (j=0; j<nbins_411; j=j+2) {
					candidate = (smspec[j] > w->min_snr) && (npk < NPK);
					if (candidate) {
						pk[npk].bin0 = j;
						pk[npk].freq0 = (j-hbins_205)*df;
						pk[npk].snr0 = 10*log10(smspec[j]) - w->snr_scaling_factor;
						npk++;
					}
				}
			} else {
				for (j=1; j<(nbins_411-1); j++) {
					candidate = (smspec[j]>smspec[j-1]) &&
								(smspec[j]>smspec[j+1]) &&
								(npk<NPK);
					if (candidate) {
						pk[npk].bin0 = j;
						pk[npk].freq0 = (j-hbins_205)*df;
						pk[npk].snr0 = 10*log10(smspec[j]) - w->snr_scaling_factor;
						npk++;
					}
				}
			}
			//wdprintf("initial npk %d/%d\n", npk, NPK);
			TRY_YIELD_DECODE;
	
			// Don't waste time on signals outside of the range [fmin,fmax].
			i=0;
			for (pki=0; pki < npk; pki++) {
				if (pk[pki].freq0 >= FMIN && pk[pki].freq0 <= FMAX) {
					pk[i].ignore = false;
					pk[i].bin0 = pk[pki].bin0;
					pk[i].freq0 = pk[pki].freq0;
					pk[i].snr0 = pk[pki].snr0;
					i++;
				}
			}
			npk = i;
			//wdprintf("freq range limited npk %d\n", npk);
			TRY_YIELD_DECODE;
			
			// only look at a limited number of strong peaks
			qsort(pk, npk, sizeof(pk_t), snr_comp);		// sort in decreasing snr order
			if (npk > MAX_NPK) npk = MAX_NPK;
			//wdprintf("MAX_NPK limited npk %d/%d\n", npk, MAX_NPK);
		
			// remember our frequency-sorted index
			qsort(pk, npk, sizeof(pk_t), freq_comp);
			for (pki=0; pki < npk; pki++) {
				pk[pki].freq_idx = pki;
			}
		
			// send peak info to client in increasing frequency order
			memcpy(pk_freq, pk, sizeof(pk));
		
			for (pki=0; pki < npk; pki++) {
				sprintf(pk_freq[pki].snr_call, "%d", (int) roundf(pk_freq[pki].snr0));
			}
		
			wspr_send_peaks(w, pk_freq, npk);
			
			// keep 'pk' in snr order so strong sigs are processed first in case we run out of decode time
			qsort(pk, npk, sizeof(pk_t), snr_comp);		// sort in decreasing snr order
		}
        
        int valid_peaks = 0;
        for (pki=0; pki < npk; pki++) {
        	if (pk[pki].ignore) continue;
        	valid_peaks++;
        }
        wdprintf("PASS %d npeaks=%d valid_peaks=%d maxcycles=%d jig_range=%+d..%d jig_step=%d ---------------------------------------------\n",
        	ipass+1, npk, valid_peaks, maxcycles, jig_range/2, -jig_range/2, iifac);

        /* Make coarse estimates of shift (DT), freq, and drift
         
         * Look for time offsets up to +/- 8 symbols (about +/- 5.4 s) relative
         to nominal start time, which is 2 seconds into the file
         
         * Calculates shift relative to the beginning of the file
         
         * Negative shifts mean that signal started before start of file
         
         * The program prints DT = shift-2 s
         
         * Shifts that cause sync vector to fall off of either end of the data
         vector are accommodated by "partial decoding", such that missing
         symbols produce a soft-decision symbol value of 128
         
         * The frequency drift model is linear, deviation of +/- drift/2 over the
         span of 162 symbols, with deviation equal to 0 at the center of the
         signal vector.
         */

        int idrift,ifr,if0,ifd,k0;
        int kindex;
        float smax,ss,power,p0,p1,p2,p3;

        for (pki=0; pki < npk; pki++) {			//For each candidate...
        	pk_t *p = &pk[pki];
            smax = -1e30;
            if0 = p->freq0/df+SPS;

			#if defined(MORE_EFFORT)
				if (ipass != 0 && p->ignore)
					continue;
			#endif
			
            for (ifr=if0-2; ifr<=if0+2; ifr++) {                      //Freq search
                for( k0=-10; k0<22; k0++) {                             //Time search
                    for (idrift=-maxdrift; idrift<=maxdrift; idrift++) {  //Drift search
                        ss=0.0;
                        power=0.0;
                        for (k=0; k<NSYM_162; k++) {				//Sum over symbols
                            ifd=ifr+((float)k-FHSYM_81)/FHSYM_81*( (float)idrift )/(2.0*df);
                            kindex=k0+2*k;
                            if( kindex < nffts ) {
								p0=w->pwr_samp[w->decode_ping_pong][ifd-3][kindex];
								p1=w->pwr_samp[w->decode_ping_pong][ifd-1][kindex];
								p2=w->pwr_samp[w->decode_ping_pong][ifd+1][kindex];
								p3=w->pwr_samp[w->decode_ping_pong][ifd+3][kindex];
                                
                                p0=sqrt(p0);
                                p1=sqrt(p1);
                                p2=sqrt(p2);
                                p3=sqrt(p3);
                                
                                ss=ss+(2*pr3[k]-1)*((p1+p3)-(p0+p2));
                                power=power+p0+p1+p2+p3;
                            }
                        }
                        sync1=ss/power;
                        if( sync1 > smax ) {                  //Save coarse parameters
                            smax=sync1;
							p->shift0=HSPS*(k0+1);
							p->drift0=idrift;
							p->freq0=(ifr-SPS)*df;
							p->sync0=sync1;
                        }
						//wdprintf("drift %d  k0 %d  sync %f\n",idrift,k0,smax);
                    }
					TRY_YIELD_DECODE;
                }
            }
			wdprintf("npeak     #%02d %6.1f snr  %9.6f (%7.2f) freq  %4.1f drift  %5d shift  %6.3f sync  %3d bin\n",
				pki, p->snr0, w->dialfreq_MHz+(bfo+p->freq0)/1e6, w->cf_offset+p->freq0, p->drift0, p->shift0, p->sync0, p->bin0);
        }

        /*
         Refine the estimates of freq, shift using sync as a metric.
         Sync is calculated such that it is a float taking values in the range
         [0.0,1.0].
         
         Function sync_and_demodulate has three modes of operation
         mode is the last argument:
         
         0 = no frequency or drift search. find best time lag.
         1 = no time lag or drift search. find best frequency.
         2 = no frequency or time lag search. Calculate soft-decision
         symbols using passed frequency and shift.
         
         NB: best possibility for OpenMP may be here: several worker threads
         could each work on one candidate at a time.
         */
		int candidates = 0;
        ndecodes_pass = 0;
		
        for (pki=0; pki < npk && !w->abort_decode; pki++) {
        	bool f_decoded = false, f_delete = false, f_image = false, f_decoding = false;

        	pk_t *p = &pk[pki];
			u4_t decode_start = timer_ms();

			#if defined(MORE_EFFORT)
				if (ipass != 0 && p->ignore)
					continue;
			#endif
			
			candidates++;
            f1 = p->freq0;
            snr = p->snr0;
            drift1 = p->drift0;
            shift1 = p->shift0;
            sync1 = p->sync0;

			f_decoding = true;
			pk_freq[p->freq_idx].flags |= WSPR_F_DECODING;
			wspr_send_peaks(w, pk_freq, npk);
	
			wdprintf("start     #%02d %6.1f snr  %9.6f (%7.2f) freq  %4.1f drift  %5d shift  %6.3f sync\n",
				pki, snr, w->dialfreq_MHz+(bfo+f1)/1e6, w->cf_offset+f1, drift1, shift1, sync1);

            // coarse-grid lag and freq search, then if sync > minsync1 continue
            fstep=0.0; ifmin=0; ifmax=0;
            lagmin = shift1-128;
            lagmax = shift1+128;
            lagstep = 64;
            sync_and_demodulate(idat, qdat, npoints, w->symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                lagmin, lagmax, lagstep, drift1, symfac, &sync1, 0);

            fstep = 0.25; ifmin = -2; ifmax = 2;
            sync_and_demodulate(idat, qdat, npoints, w->symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                lagmin, lagmax, lagstep, drift1, symfac, &sync1, 1);

            // refine drift estimate
            fstep=0.0; ifmin=0; ifmax=0;
            float driftp,driftm,syncp,syncm;
            driftp = drift1+0.5;
            sync_and_demodulate(idat, qdat, npoints, w->symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                lagmin, lagmax, lagstep, driftp, symfac, &syncp, 1);
            
            driftm = drift1-0.5;
            sync_and_demodulate(idat, qdat, npoints, w->symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                lagmin, lagmax, lagstep, driftm, symfac, &syncm, 1);
            
            if (syncp > sync1) {
                drift1 = driftp;
                sync1 = syncp;
            } else
            
            if (syncm > sync1) {
                drift1 = driftm;
                sync1 = syncm;
            }

			wdprintf("coarse    #%02d %6.1f snr  %9.6f (%7.2f) freq  %4.1f drift  %5d shift  %6.3f sync\n",
				pki, snr, w->dialfreq_MHz+(bfo+f1)/1e6, w->cf_offset+f1, drift1, shift1, sync1);

            // fine-grid lag and freq search
			bool r_minsync1 = (sync1 > minsync1);

            if (r_minsync1) {
                lagmin = shift1-32; lagmax = shift1+32; lagstep = 16;
                sync_and_demodulate(idat, qdat, npoints, w->symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                    lagmin, lagmax, lagstep, drift1, symfac, &sync1, 0);
            
                // fine search over frequency
                fstep = 0.05; ifmin = -2; ifmax = 2;
                sync_and_demodulate(idat, qdat, npoints, w->symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                lagmin, lagmax, lagstep, drift1, symfac, &sync1, 1);
            } else {
            	p->ignore = true;
				wdprintf("MINSYNC1  #%02d\n", pki);
				f_delete = true;
            }
            
            int idt=0, ii=0, jiggered_shift;
            float y, sq, rms;
        	int r_decoded = 0;
        	bool r_tooWeak = true;
            
            // ii: 0 +1 -1 +2 -2 +3 -3 ... (*iifac)
            // ii always covers jig_range, stepped by iifac resolution
            while (!w->abort_decode && r_minsync1 && !r_decoded && idt <= (jig_range/iifac)) {
                ii = (idt+1)/2;
                if ((idt&1) == 1) ii = -ii;
                ii = iifac*ii;
                jiggered_shift = shift1+ii;
                
                // Use mode 2 to get soft-decision symbols
                sync_and_demodulate(idat, qdat, npoints, w->symbols, &f1, ifmin, ifmax, fstep,
                                    &jiggered_shift, lagmin, lagmax, lagstep, drift1, symfac,
                                    &sync1, 2);

                sq = 0.0;
                for (i=0; i<NSYM_162; i++) {
                    y = (float) w->symbols[i] - 128.0;
                    sq += y*y;
                }
                rms = sqrt(sq/FNSYM_162);

				bool weak = true;
                if ((sync1 > minsync2) && (rms > minrms)) {
                    deinterleave(w->symbols);
                    
                    if (w->stackdecoder) {
                        r_decoded = jelinek(&metric, &cycles, w->decdata, w->symbols, NBITS,
											WSPR_STACKSIZE, w->stack, mettab, maxcycles);
                    } else {
                        r_decoded = fano(&metric, &cycles, &maxnp, w->decdata, w->symbols, NBITS,
										mettab, delta, maxcycles);
                    }

                    r_tooWeak = weak = false;
                }
                
				wdprintf("jig <>%3d #%02d %6.1f snr  %9.6f (%7.2f) freq  %4.1f drift  %5d(%+4d) shift  %6.3f sync  %4.1f rms",
					idt, pki, snr, w->dialfreq_MHz+(bfo+f1)/1e6, w->cf_offset+f1, drift1, jiggered_shift, ii, sync1, rms);
				if (!weak) {
					wprintf("  %4u metric  %3u cycles\n", metric, cycles);
				} else {
					if (sync1 <= minsync2) wprintf("  SYNC-WEAK");
					if (rms <= minrms) wprintf("  RMS-WEAK");
					wprintf("\n");
				}
				
                idt++;
                if (w->quickmode) break;
            }
            
            bool r_timeUp = (!w->abort_decode && r_minsync1 && !r_decoded && idt > (jig_range/iifac));
            int r_valid = 0;
            
            //if (r_timeUp && r_tooWeak && iifac == 1) {
            if (r_timeUp && r_tooWeak) {
				wdprintf("NO CHANGE #%02d\n", pki);
				p->ignore = true;	// situation not going to get any better
				f_delete = true;
			}
			
			// result priority: abort, minSync1, tooWeak, noChange, timeUp, image, decoded
            
            if (r_decoded) {
                ndecodes_pass++;
                p->ignore = true;
                
                // Unpack the decoded message, update the hashtable, apply
                // sanity checks on grid and power, and return
                // call_loc_pow string and also callsign (for de-duping).
                int dBm;
                r_valid = unpk_(w->decdata, w->call_loc_pow, w->callsign, w->grid, &dBm);

                // subtract even on last pass
                #ifdef SUBTRACT_SIGNAL
					if (w->subtraction && (ipass < npasses) && r_valid > 0) {
						if (get_wspr_channel_symbols(w->call_loc_pow, w->channel_symbols)) {
							subtract_signal2(idat, qdat, npoints, f1, shift1, drift1, w->channel_symbols);
						} else {
							break;
						}
						
					}
                #endif

                // Remove dupes (same callsign and freq within 3 Hz) and images
                bool r_dupe = false;
                if (r_valid > 0) {
                	bool hash = (strcmp(w->callsign, "...") == 0);
                	for (i=0; i < uniques; i++) {
                		decode_t *dp = &w->deco[i];
						bool match = (strcmp(w->callsign, dp->call) == 0);
						bool close = (fabs(f1 - dp->freq) < 3.0);
						if ((match && close) || (match && !hash)) {
							if (close) {
								f_delete = true;
							} else {
								f_image = true;
							}
							wdprintf("%s     #%02d  with #%02d %s, %.3f secs\n", f_image? "IMAGE" : "DUPE ",
								pki, i, dp->call, (float)(timer_ms()-decode_start)/1e3);
							r_dupe = true;
							break;
						}
					}
				}

				if (r_valid <= 0 && !r_dupe) {
					if (r_valid < 0) {
						wdprintf("UNPK ERR! #%02d  error code %d, %.3f secs\n",
							pki, r_valid, (float)(timer_ms()-decode_start)/1e3);
					} else {
						wdprintf("NOT VALID #%02d  %.3f secs\n",
							pki, (float)(timer_ms()-decode_start)/1e3);
					}
					f_delete = true;
				}

                if (r_valid > 0 && !r_dupe) {
                	f_decoded = true;

                    if (w->wspr_type == WSPR_TYPE_15MIN) {
                        freq_print = w->dialfreq_MHz + (bfo+112.5+f1/8.0)/1e6;
                        dt_print = shift1*8*dt-1.0;
                    } else {
                        freq_print = w->dialfreq_MHz + (bfo+f1)/1e6;
                        dt_print = shift1*dt-1.0;
                    }

					struct tm tm;
					gmtime_r(&w->utc[w->decode_ping_pong], &tm);
            
					wdprintf("TYPE%d %02d%02d %3.0f %4.1f %10.6f %2d %-s %4s %2d [%s] in %.3f secs --------------------------------------------------------------------\n",
					   r_valid, tm.tm_hour, tm.tm_min, snr, dt_print, freq_print, (int) drift1,
					   w->callsign, w->grid, dBm, w->call_loc_pow, (float)(timer_ms()-decode_start)/1e3);
					
					double watts, factor;
					char *W_s;
					const char *units;
	
					watts = pow(10.0, (dBm - 30)/10.0);
					if (watts >= 1.0) {
						factor = 1.0;
						units = "W";
					} else
					if (watts >= 1e-3) {
						factor = 1e3;
						units = "mW";
					} else
					if (watts >= 1e-6) {
						factor = 1e6;
						units = "uW";
					} else
					if (watts >= 1e-9) {
						factor = 1e9;
						units = "nW";
					} else {
						factor = 1e12;
						units = "pW";
					}
					
					watts *= factor;
					if (watts < 10.0)
						asprintf(&W_s, "%.1f %s", watts, units);
					else
						asprintf(&W_s, "%.0f %s", watts, units);
					
					// 			Call   Grid    km  dBm
					// TYPE1	c6cccc g4gg kkkkk  ppp (s)
					// TYPE2	...                ppp (s)		; no hash yet
					// TYPE2	ppp/c6cccc         ppp (s)		; 1-3 char prefix
					// TYPE2	c6cccc/ss          ppp (s)		; 1-2 char suffix
					// TYPE3	...  g6gggg kkkkk  ppp (s)		; no hash yet
					// TYPE3	ppp/c6cccc g6gggg kkkkk ppp (s)	; worst case, just let the fields float

					if (r_valid == 1)		// TYPE1
						fprintf(w->demod_pipe,
							"%02d%02d %3.0f %4.1f %9.6f %2d  "
							"<a href='https://www.qrz.com/lookup/%s' target='_blank'>%-6s</a> "
							"<a href='http://www.levinecentral.com/ham/grid_square.php?Grid=%s' target='_blank'>%s</a> "
							"%5d  %3d (%s)\n",
							tm.tm_hour, tm.tm_min, snr, dt_print, freq_print, (int) drift1,
							w->callsign, w->callsign, w->grid, w->grid,
							(int) grid_to_distance_km(w->grid), dBm, W_s);
					else
					
					if (r_valid == 2) {	// TYPE2
						if (strcmp(w->callsign, "...") == 0)
						    fprintf(w->demod_pipe,
								"%02d%02d %3.0f %4.1f %9.6f %2d  ...                %3d (%s)\n",
								tm.tm_hour, tm.tm_min, snr, dt_print, freq_print, (int) drift1, dBm, W_s);
						else
						    fprintf(w->demod_pipe,
								"%02d%02d %3.0f %4.1f %9.6f %2d  "
								"<a href='https://www.qrz.com/lookup/%s' target='_blank'>%-17s</a>"
								"  %3d (%s)\n",
							   //kkkkk--ppp
								tm.tm_hour, tm.tm_min, snr, dt_print, freq_print, (int) drift1,
								w->callsign, w->callsign, dBm, W_s);
					} else
					
					if (r_valid == 3) {	// TYPE3
						if (strcmp(w->callsign, "...") == 0)
						    fprintf(w->demod_pipe,
								"%02d%02d %3.0f %4.1f %9.6f %2d  "
								"...  "
								"<a href='http://www.levinecentral.com/ham/grid_square.php?Grid=%s' target='_blank'>%6s</a> "
								"%5d  %3d (%s)\n",
								tm.tm_hour, tm.tm_min, snr, dt_print, freq_print, (int) drift1,
								w->grid, w->grid, (int) grid_to_distance_km(w->grid), dBm, W_s);
						else
						    fprintf(w->demod_pipe,
								"%02d%02d %3.0f %4.1f %9.6f %2d  "
								"<a href='https://www.qrz.com/lookup/%s' target='_blank'>%s</a> "
								"<a href='http://www.levinecentral.com/ham/grid_square.php?Grid=%s' target='_blank'>%s</a> "
								"%d %d (%s)\n",
								tm.tm_hour, tm.tm_min, snr, dt_print, freq_print, (int) drift1,
								w->callsign, w->callsign, w->grid, w->grid,
								(int) grid_to_distance_km(w->grid), dBm, W_s);
					}
					
					free(W_s);
					fflush(w->demod_pipe);
					
					decode_t *dp = &w->deco[uniques];
                    strcpy(dp->call, w->callsign);
                    dp->freq = f1;
                    dp->hour = tm.tm_hour;
                    dp->min = tm.tm_min;
                    dp->snr = snr;
                    dp->dt_print = dt_print;
                    dp->freq_print = freq_print;
                    dp->drift1 = drift1;
                    strcpy(dp->c_l_p, w->call_loc_pow);
                    uniques++;
				} else {
					r_valid = 0;
                }
			} else {
				if (r_timeUp) {
					wdprintf("TIME UP   #%02d %.3f secs\n", pki, (float)(timer_ms()-decode_start)/1e3);
					#if defined(MORE_EFFORT)
						f_decoding = false;
					#else
						f_delete = true;
					#endif
				}
            }	// decoded

			if (f_decoded) {
				strcpy(pk_freq[p->freq_idx].snr_call, w->callsign);
				pk_freq[p->freq_idx].flags |= WSPR_F_DECODED;
			} else
			if (f_delete) {
				pk_freq[p->freq_idx].flags |= WSPR_F_DELETE;
			} else
			if (f_image) {
				strcpy(pk_freq[p->freq_idx].snr_call, "image");
				pk_freq[p->freq_idx].flags |= WSPR_F_IMAGE;
			} else
			if (!f_decoding) {
				pk_freq[p->freq_idx].flags &= ~WSPR_F_DECODING;
			}
			
			wspr_send_peaks(w, pk_freq, npk);
        }	// peak list
        
		if (candidates == 0)
			break;		// nothing left to do
			
    }	// passes

	// when finished delete any unresolved peaks
	for (i=0; i < npk; i++) {
		pk_t *p = &pk[i];
		if (!(pk_freq[p->freq_idx].flags & WSPR_F_DECODED))
			pk_freq[p->freq_idx].flags |= WSPR_F_DELETE;
	}
	wspr_send_peaks(w, pk_freq, npk);

	// upload spots at the end of the decoding when there is less load on wsprnet.org
	// FIXME OpenWebRX doesn't upload yet
	for (i = 0; i < uniques; i++) {
		decode_t *dp = &w->deco[i];
		ext_send_msg_encoded(w->rx_chan, WSPR_DEBUG_MSG, "EXT", "WSPR_UPLOAD",
			"%02d%02d %3.0f %4.1f %9.6f %2d %s",
			dp->hour, dp->min, dp->snr, dp->dt_print, dp->freq_print, (int) dp->drift1, dp->c_l_p);
		wdprintf("WSPR_UPLOAD U%d/%d "
			"%02d%02d %3.0f %4.1f %9.6f %2d %s" "\n", i, uniques,
			dp->hour, dp->min, dp->snr, dp->dt_print, dp->freq_print, (int) dp->drift1, dp->c_l_p);
	}
}
