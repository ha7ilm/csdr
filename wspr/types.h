#pragma once

#include <stdlib.h>

typedef unsigned long long	u64_t;
typedef unsigned int        u4_t;
typedef unsigned char       u1_t;
typedef unsigned short      u2_t;

typedef signed long long	s64_t;
typedef signed int			s4_t;
typedef signed short        s2_t;
typedef signed char         s1_t;

typedef void (*func_t)();
typedef void (*funcP_t)(void *);
typedef int (*funcPR_t)(void *);

#define TO_VOID_PARAM(p)    ((void *) (long) (p))
#define FROM_VOID_PARAM(p)  ((long) (p))

#define U1(v) ((u1_t) (v))
#define S1(v) ((s1_t) (v))
#define U2(v) ((u2_t) (v))
#define S2(v) ((s2_t) (v))
#define U4(v) ((u4_t) (v))
#define S4(v) ((s4_t) (v))
#define U8(v) ((u64_t) (v))
#define S8(v) ((s64_t) (v))

#define S16x4_S64(a,b,c,d)	S8( (U8(a)<<48) | (U8(b)<<32) | (U8(c)<<16) | U8(d) )
#define S14_16(w)			S2( U2(w) | ( ((U2(w)) & 0x2000)? 0xc000:0 ) )
#define S14_32(w)			S4( S2( U2(w) | ((U2(w) & 0x2000)? 0xffffc000:0 ) ) )
#define S24_8_16(h8,l16)	S4( (U1(h8)<<16) | U2(l16) | ((U1(h8) & 0x80)? 0xff000000:0) )
#define S24_16_8(h16,l8)	S4( (U2(h16)<<8) | U1(l8) | ((U2(h16) & 0x8000)? 0xff000000:0) )

#define B2I(bytes)			(((bytes)+3)/4)
#define I2B(ints)			((ints)*4)
#define B2S(bytes)			(((bytes)+1)/2)
#define S2B(shorts)			((shorts)*2)

#define	B3(i)				(((i) >> 24) & 0xff)
#define	B2(i)				(((i) >> 16) & 0xff)
#define	B1(i)				(((i) >>  8) & 0xff)
#define	B0(i)				(((i) >>  0) & 0xff)

#define	FLIP32(i)			((B0(i) << 24) | (B1(i) << 16) | (B2(i) << 8) | (B3(i) << 0))
#define	FLIP16(i)			((B0(i) << 8) | (B1(i) << 0))

#ifndef TRUE
 #define TRUE 1
#endif

#ifndef FALSE
 #define FALSE 0
#endif

#ifndef __cplusplus
 typedef	unsigned char	bool;
 #define true TRUE
 #define false FALSE
#endif

#define	NOT_FOUND	-1

#define	ARRAY_LEN(x)	((int) (sizeof (x) / sizeof ((x) [0])))

#define	K		1024
#define	M		(K*K)
#define	B		(M*K)

#define	MHz		1000000
#define	kHz		1000

#define MAX(a,b) ((a)>(b)?(a):(b))
#define max(a,b) MAX(a,b)
#define MIN(a,b) ((a)<(b)?(a):(b))
#define min(a,b) MIN(a,b)

#define SI_CLAMP(a,n) ( ((a) > ((n)-1))? ((n)-1) : ( ((a) < -(n))? -(n) : (a) ) )

#define	STRINGIFY(x) #x
#define	STRINGIFY_DEFINE(x) STRINGIFY(x)	// indirection needed for a -Dx=y define
#define	CAT_STRING(x,y) x y			// just a reminder of how this is done: "foo" "bar"
#define	CAT_DEFINE_VAR(x,y) x ## y	// just a reminder of how this is done: foo ## bar

// documentation assistance
#define SPACE_FOR_NULL 1
