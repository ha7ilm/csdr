/*
 Soft-decision stack-based sequential decoder for K=32 r=1/2
 convolutional code. This code implements the "stack-bucket" algorithm
 described in:
 "Fast Sequential Decoding Algorithm Using a Stack", F. Jelinek
 
 The ENCODE macro from Phil Karn's (KA9Q) Fano decoder is used.
 
 Written by Steve Franke, K9AN for WSJT-X (July 2015)
 */

#include "jelinek.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> /* memset */

#include "fano.h"

/* WSPR uses the Layland-Lushbaugh code
 * Nonsystematic, non-quick look-in, dmin=?, dfree=?
 */
#define	POLY1	0xf2d05351
#define	POLY2	0xe4613c47

//Decoder - returns 1 on success, 0 on timeout
int jelinek(
            unsigned int *metric,	/* Final path metric (returned value) */
            unsigned int *cycles,	/* Cycle count (returned value) */
            unsigned char *data,	/* Decoded output data */
            unsigned char *symbols,	/* Raw deinterleaved input symbols */
            unsigned int nbits,	/* Number of output bits */
            unsigned int stacksize,
            struct snode *stack,
            int mettab[2][256],	/* Metric table, [sent sym][rx symbol] */
            unsigned int maxcycles)/* Decoding timeout in cycles per bit */
{
    
    // Compute branch metrics for each symbol pair
    // The sequential decoding algorithm only uses the metrics, not the
    // symbol values.
    unsigned int i;
    long int metrics[81][4];
    for(i=0; i<nbits; i++){
        metrics[i][0] = mettab[0][symbols[0]] + mettab[0][symbols[1]];
        metrics[i][1] = mettab[0][symbols[0]] + mettab[1][symbols[1]];
        metrics[i][2] = mettab[1][symbols[0]] + mettab[0][symbols[1]];
        metrics[i][3] = mettab[1][symbols[0]] + mettab[1][symbols[1]];
        symbols += 2;
    }
    
    // zero the stack
    // jksd very expensive given large enough stacksize
    memset(stack,0,stacksize*sizeof(struct snode));
    
    // initialize the loop variables
    unsigned int lsym, ntail=31;
    uint64_t encstate=0;
    unsigned int nbuckets=1000;
    unsigned int low_bucket=nbuckets-1; //will be set on first run-through
    unsigned int high_bucket=0;
    unsigned int *buckets, bucket;
    buckets = (unsigned int *) malloc(nbuckets*sizeof(unsigned int));
    memset(buckets,0,nbuckets*sizeof(unsigned int));
    unsigned int ptr=1;
    unsigned int stackptr=1; //pointer values of 0 are reserved (they mean that a bucket is empty)
    unsigned int depth=0, nbits_minus_ntail=nbits-ntail;
    unsigned int stacksize_minus_1=stacksize-1;
    long int totmet0, totmet1, gamma=0;
    
    unsigned int ncycles=maxcycles*nbits;
    /********************* Start the stack decoder *****************/
    for (i=1; i <= ncycles; i++) {
#ifdef DEBUG
        printf("***stackptr=%ud, depth=%d, gamma=%ld, encstate=%llx, bucket %d, low_bucket %d, high_bucket %d\n",
               stackptr, depth, gamma, encstate, bucket, low_bucket, high_bucket);
#endif
        // no need to store more than 7 bytes (56 bits) for encoder state because
        // only 50 bits are not 0's.
        if( depth < 56 ) {
            encstate=encstate<<1;
            ENCODE(lsym,encstate); // get channel symbols associated with the 0 branch
        } else {
            ENCODE(lsym,encstate<<(depth-55));
        }

        // lsym are the 0-branch channel symbols and 3^lsym are the 1-branch
        // channel symbols (due to a special property of our generator polynomials)
        totmet0 = gamma+metrics[depth][lsym];   // total metric for 0-branch daughter node
        totmet1 = gamma+metrics[depth][3^lsym]; // total metric for 1-branch daughter node
        depth++; //the depth of the daughter nodes

        bucket=(totmet0>>5)+200; //fast, but not particularly safe - totmet can be negative
        if( bucket > high_bucket ) high_bucket=bucket;
        if( bucket < low_bucket ) low_bucket=bucket;
       
        // place the 0 node on the stack, overwriting the parent (current) node
        stack[ptr].encstate=encstate;
        stack[ptr].gamma=totmet0;
        stack[ptr].depth=depth;
        stack[ptr].jpointer=buckets[bucket];
        buckets[bucket]=ptr;
        
        // if in the tail, only need to evaluate the "0" branch.
        // Otherwise, enter this "if" and place the 1 node on the stack,
        if( depth <= nbits_minus_ntail ) {
            if( stackptr < stacksize_minus_1 ) {
                stackptr++;
                ptr=stackptr;
            } else { // stack full
                while( buckets[low_bucket] == 0 ) { //write latest to where the top of the lowest bucket points
                    low_bucket++;
                }
                ptr=buckets[low_bucket];
                buckets[low_bucket]=stack[ptr].jpointer; //make bucket point to next older entry
            }

            bucket=(totmet1>>5)+200; //this may not be safe on all compilers
            if( bucket > high_bucket ) high_bucket=bucket;
            if( bucket < low_bucket ) low_bucket=bucket;
            
            stack[ptr].encstate=encstate+1;
            stack[ptr].gamma=totmet1;
            stack[ptr].depth=depth;
            stack[ptr].jpointer=buckets[bucket];
            buckets[bucket]=ptr;
        }

    // pick off the latest entry from the high bucket
        while( buckets[high_bucket] == 0 ) {
            high_bucket--;
        }
        ptr=buckets[high_bucket];
        buckets[high_bucket]=stack[ptr].jpointer;
        depth=stack[ptr].depth;
        gamma=stack[ptr].gamma;
        encstate=stack[ptr].encstate;

        // we are done if the top entry on the stack is at depth nbits
        if (depth == nbits) {
            break;
        }
    }
    
    *cycles = i+1;
    *metric =  gamma;	/* Return final path metric */

    //    printf("cycles %d stackptr=%d, depth=%d, gamma=%d, encstate=%lx\n",
    //           *cycles, stackptr, depth, *metric, encstate);
    
    for (i=0; i<7; i++) {
        data[i]=(encstate>>(48-i*8))&(0x00000000000000ff);
    }
    for (i=7; i<11; i++) {
        data[i]=0;
    }

    if(*cycles/nbits >= maxcycles) //timed out
    {
        return 0;
    }
    return 1;		//success
}
