/*

This software is part of libcsdr, a set of simple DSP routines for 
Software Defined Radio.

Copyright (c) 2015, Andras Retzler <randras@sdr.hu>
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

   Copyright 1997 Tim Kientzle.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. All advertising materials mentioning features or use of this software
   must display the following acknowledgement:
      This product includes software developed by Tim Kientzle
      and published in ``The Programmer's Guide to Sound.''
4. Neither the names of Tim Kientzle nor Addison-Wesley
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL TIM KIENTZLE OR ADDISON-WESLEY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/***********************************************************
Copyright 1992 by Stichting Mathematisch Centrum, Amsterdam, The
Netherlands.

                        All Rights Reserved

Permission to use, copy, modify, and distribute this software and its 
documentation for any purpose and without fee is hereby granted, 
provided that the above copyright notice appear in all copies and that
both that copyright notice and this permission notice appear in 
supporting documentation, and that the names of Stichting Mathematisch
Centrum or CWI not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior permission.

STICHTING MATHEMATISCH CENTRUM DISCLAIMS ALL WARRANTIES WITH REGARD TO
THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS, IN NO EVENT SHALL STICHTING MATHEMATISCH CENTRUM BE LIABLE
FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT
OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

******************************************************************/

#ifdef USE_IMA_ADPCM

#include "ima_adpcm.h"
 
const int indexAdjustTable[16] = {
   -1, -1, -1, -1,  // +0 - +3, decrease the step size
    2, 4, 6, 8,     // +4 - +7, increase the step size
   -1, -1, -1, -1,  // -0 - -3, decrease the step size
    2, 4, 6, 8,     // -4 - -7, increase the step size
};

const int _stepSizeTable[89] = {
   7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 19, 21, 23, 25, 28, 31, 34,
   37, 41, 45, 50, 55, 60, 66, 73, 80, 88, 97, 107, 118, 130, 143,
   157, 173, 190, 209, 230, 253, 279, 307, 337, 371, 408, 449, 494,
   544, 598, 658, 724, 796, 876, 963, 1060, 1166, 1282, 1411, 1552,
   1707, 1878, 2066, 2272, 2499, 2749, 3024, 3327, 3660, 4026,
   4428, 4871, 5358, 5894, 6484, 7132, 7845, 8630, 9493, 10442,
   11487, 12635, 13899, 15289, 16818, 18500, 20350, 22385, 24623,
   27086, 29794, 32767
};

static inline short ImaAdpcmDecode(unsigned char deltaCode, ima_adpcm_state_t* state) {
   // Get the current step size
   int step = _stepSizeTable[state->index];

   // Construct the difference by scaling the current step size
   // This is approximately: difference = (deltaCode+.5)*step/4
   int difference = step>>3;
   if ( deltaCode & 1 ) difference += step>>2;
   if ( deltaCode & 2 ) difference += step>>1;
   if ( deltaCode & 4 ) difference += step;
   if ( deltaCode & 8 ) difference = -difference;

   // Build the new sample
   state->previousValue += difference;
   if (state->previousValue > 32767) state->previousValue = 32767;
   else if (state->previousValue < -32768) state->previousValue = -32768;

   // Update the step for the next sample
   state->index += indexAdjustTable[deltaCode];
   if (state->index < 0) state->index = 0;
   else if (state->index > 88) state->index = 88;

   return state->previousValue;
}

static inline unsigned char ImaAdpcmEncode(short sample, ima_adpcm_state_t* state) {
   int diff = sample - state->previousValue;
   int step = _stepSizeTable[state->index];
   int deltaCode = 0;

   // Set sign bit
   if (diff < 0) { deltaCode = 8; diff = -diff; }

   // This is essentially deltaCode = (diff<<2)/step,
   // except the roundoff is handled differently.
   if ( diff >= step ) {  deltaCode |= 4;  diff -= step;  }
   step >>= 1;
   if ( diff >= step ) {  deltaCode |= 2;  diff -= step;  }
   step >>= 1;
   if ( diff >= step ) {  deltaCode |= 1;  diff -= step;  }

   ImaAdpcmDecode(deltaCode,state);  // update state
   return deltaCode;
}

ima_adpcm_state_t encode_ima_adpcm_i16_u8(short* input, unsigned char* output, int input_length, ima_adpcm_state_t state)
{
	int k=0;
	for(int i=0;i<input_length/2;i++) 
	{
		output[k]=ImaAdpcmEncode(input[2*i],&state);
		output[k++]|=ImaAdpcmEncode(input[2*i+1],&state)<<4;
	}
	return state;
}

ima_adpcm_state_t decode_ima_adpcm_u8_i16(unsigned char* input, short* output, int input_length, ima_adpcm_state_t state)
{
	int k=0;
	for(int i=0;i<input_length;i++) 
	{
		output[k++]=ImaAdpcmDecode(input[i]&0xf,&state);
		output[k++]=ImaAdpcmDecode( (input[i]>>4)&0xf,&state);
	}
	return state;
}

#endif
