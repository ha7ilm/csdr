/*
 This file is part of wsprd.
 
 File name: fano.h

 Description: Header file for sequential Fano decoder.

 Copyright 1994, Phil Karn, KA9Q
 Minor modifications by Joe Taylor, K1JT
*/

#ifndef FANO_H
#define FANO_H

int fano(unsigned int *metric, unsigned int *cycles, unsigned int *maxnp,
	unsigned char *data,unsigned char *symbols, unsigned int nbits,
	 int mettab[2][256],int delta,unsigned int maxcycles);

int encode(unsigned char *symbols,unsigned char *data,unsigned int nbytes);

extern unsigned char Partab[];

/* Convolutional encoder macro. Takes the encoder state, generates
 * a rate 1/2 symbol pair and stores it in 'sym'. The symbol generated from
 * POLY1 goes into the 2-bit of sym, and the symbol generated from POLY2
 * goes into the 1-bit.
 */
#define	ENCODE(sym,encstate){\
unsigned long _tmp;\
\
_tmp = (encstate) & POLY1;\
_tmp ^= _tmp >> 16;\
(sym) = Partab[(_tmp ^ (_tmp >> 8)) & 0xff] << 1;\
_tmp = (encstate) & POLY2;\
_tmp ^= _tmp >> 16;\
(sym) |= Partab[(_tmp ^ (_tmp >> 8)) & 0xff];\
}

#endif
