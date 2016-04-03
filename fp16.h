/*
This software is part of libcsdr, a set of simple DSP routines for 
Software Defined Radio.

Copyright (c) 2016, Andras Retzler <randras@sdr.hu>
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

// This code originates from: https://gist.githubusercontent.com/rygorous/2156668/raw/ef8408efac2ff0db549252883dd4c99dddfcc929/gistfile1.cpp
// It is the great work of Fabian Giesen.

// float->half variants.
// by Fabian "ryg" Giesen.
//
// I hereby place this code in the public domain, as per the terms of the
// CC0 license:
//
//   https://creativecommons.org/publicdomain/zero/1.0/
//
// float_to_half_full: This is basically the ISPC stdlib code, except
// I preserve the sign of NaNs (any good reason not to?)
//
// float_to_half_fast: Almost the same, with some unnecessary cases cut.
//
// float_to_half_fast2: This is where it gets a bit weird. See lengthy
// commentary inside the function code. I'm a bit on the fence about two
// things:
// 1. This *will* behave differently based on whether flush-to-zero is
//    enabled or not. Is this acceptable for ISPC?
// 2. I'm a bit on the fence about NaNs. For half->float, I opted to extend
//    the mantissa (preserving both qNaN and sNaN contents) instead of always
//    returning a qNaN like the original ISPC stdlib code did. For float->half
//    the "natural" thing would be just taking the top mantissa bits, except
//    that doesn't work; if they're all zero, we might turn a sNaN into an
//    Infinity (seriously bad!). I could test for this case and do a sticky
//    bit-like mechanism, but that's pretty ugly. Instead I go with ISPC
//    std lib behavior in this case and just return a qNaN - not quite symmetric
//    but at least it's always safe. Any opinions?
//
// I'll just go ahead and give "fast2" the same treatment as my half->float code,
// but if there's concerns with the way it works I might revise it later, so watch
// this spot.
//
// float_to_half_fast3: Bitfields removed. Ready for SSE2-ification :)
//
// float_to_half_SSE2: Exactly what it says on the tin. Beware, this works slightly
// differently from float_to_half_fast3 - the clamp and bias steps in the "normal" path
// are interchanged, since I get "minps" on every SSE2 target, but "pminsd" only for
// SSE4.1 targets. This code does what it should and is remarkably short, but the way
// it ended up working is "nonobvious" to phrase it politely.
//
// approx_float_to_half: Simpler (but less accurate) version that matches the Fox
// toolkit float->half conversions: http://blog.fox-toolkit.org/?p=40 - note that this
// also (incorrectly) translates some sNaNs into infinity, so be careful!
//
// approx_float_to_half_SSE2: SSE2 version of above.
//
// ----
//
// UPDATE 2016-01-25: Now also with a variant that implements proper round-to-nearest-even.
// It's a bit more expensive and has seen less tweaking than the other variants. On the
// plus side, it doesn't produce subnormal FP32 values as part of generating subnormal
// FP16 values, so the performance is a lot more consistent.
//
// float_to_half_rtne_full: Unoptimized round-to-nearest-break-ties-to-even reference
// implementation.
// 
// float_to_half_fast3_rtne: Variant of float_to_half_fast3 that performs round-to-
// nearest-even.
//
// float_to_half_rtne_SSE2: SSE2 implementation of float_to_half_fast3_rtne.
//
// All three functions have been exhaustively tested to produce the same results on
// all 32-bit floating-point numbers with SSE arithmetic in round-to-nearest-even mode.
// No guarantees for what happens with other rounding modes! (See testbed code.)
//
// ----
//
// Oh, and enumerating+testing all 32-bit floats takes some time, especially since
// we will snap a significant fraction of the overall FP32 range to denormals, not
// exactly a fast operation. There's a reason this one prints regular progress
// reports. You've been warned.


#include <stdio.h>
#include <stdlib.h>
#include <emmintrin.h>

typedef unsigned int uint;

union FP32_u
{
    uint u;
    float f;
    struct
    {
        uint Mantissa : 23;
        uint Exponent : 8;
        uint Sign : 1;
    };
};

union FP16_u
{
    unsigned short u;
    struct
    {
        uint Mantissa : 10;
        uint Exponent : 5;
        uint Sign : 1;
    };
};

typedef union FP32_u FP32;
typedef union FP16_u FP16;


FP16 float_to_half_full(FP32 f);
FP16 float_to_half_full_rtne(FP32 f);
FP16 float_to_half_fast(FP32 f);
FP16 float_to_half_fast2(FP32 f);
FP16 float_to_half_fast3(FP32 f);
FP16 float_to_half_fast3_rtne(FP32 f);
FP16 approx_float_to_half(FP32 f);
__m128i float_to_half_SSE2(__m128 f);
__m128i float_to_half_rtne_SSE2(__m128 f);
__m128i approx_float_to_half_SSE2(__m128 f);
void f16_generatetables();
uint float_to_half_foxtk(uint f);
FP32 half_to_float(FP16 h);
FP32 half_to_float_lit(unsigned short u);

void convert_f_f16(float* input, short* output, int input_size);
void convert_f16_f(short* input, float* output, int input_size);

