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

#include "fp16.h"

void convert_f_f16(float* input, short* output, int input_size)
{
	FP32 f32;
	for(int i=0;i<input_size;i++) //@convert_f_f16
	{
		f32.f=input[i];
		output[i]=float_to_half_full(f32).u; 
	}
}

void convert_f16_f(short* input, float* output, int input_size)
{
	FP16 f16;
	for(int i=0;i<input_size;i++) //@convert_f16_f
	{ 
		f16.u=input[i];
		output[i]=half_to_float(f16).f; 
	}
}

// Original ISPC reference version; this always rounds ties up.
FP16 float_to_half_full(FP32 f)
{
    FP16 o = { 0 };

    // Based on ISPC reference code (with minor modifications)
    if (f.Exponent == 0) // Signed zero/denormal (which will underflow)
        o.Exponent = 0;
    else if (f.Exponent == 255) // Inf or NaN (all exponent bits set)
    {
        o.Exponent = 31;
        o.Mantissa = f.Mantissa ? 0x200 : 0; // NaN->qNaN and Inf->Inf
    }
    else // Normalized number
    {
        // Exponent unbias the single, then bias the halfp
        int newexp = f.Exponent - 127 + 15;
        if (newexp >= 31) // Overflow, return signed infinity
            o.Exponent = 31;
        else if (newexp <= 0) // Underflow
        {
            if ((14 - newexp) <= 24) // Mantissa might be non-zero
            {
                uint mant = f.Mantissa | 0x800000; // Hidden 1 bit
                o.Mantissa = mant >> (14 - newexp);
                if ((mant >> (13 - newexp)) & 1) // Check for rounding
                    o.u++; // Round, might overflow into exp bit, but this is OK
            }
        }
        else
        {
            o.Exponent = newexp;
            o.Mantissa = f.Mantissa >> 13;
            if (f.Mantissa & 0x1000) // Check for rounding
                o.u++; // Round, might overflow to inf, this is OK
        }
    }

    o.Sign = f.Sign;
    return o;
}

// Same as above, but with full round-to-nearest-even.
FP16 float_to_half_full_rtne(FP32 f)
{
    FP16 o = { 0 };

    // Based on ISPC reference code (with minor modifications)
    if (f.Exponent == 0) // Signed zero/denormal (which will underflow)
        o.Exponent = 0;
    else if (f.Exponent == 255) // Inf or NaN (all exponent bits set)
    {
        o.Exponent = 31;
        o.Mantissa = f.Mantissa ? 0x200 : 0; // NaN->qNaN and Inf->Inf
    }
    else // Normalized number
    {
        // Exponent unbias the single, then bias the halfp
        int newexp = f.Exponent - 127 + 15;
        if (newexp >= 31) // Overflow, return signed infinity
            o.Exponent = 31;
        else if (newexp <= 0) // Underflow
        {
            if ((14 - newexp) <= 24) // Mantissa might be non-zero
            {
                uint mant = f.Mantissa | 0x800000; // Hidden 1 bit
                uint shift = 14 - newexp;
                o.Mantissa = mant >> shift;

                uint lowmant = mant & ((1 << shift) - 1);
                uint halfway = 1 << (shift - 1);

                if (lowmant >= halfway) // Check for rounding
                {
                    if (lowmant > halfway || (o.Mantissa & 1)) // if above halfway point or unrounded result is odd
                        o.u++; // Round, might overflow into exp bit, but this is OK
                }
            }
        }
        else
        {
            o.Exponent = newexp;
            o.Mantissa = f.Mantissa >> 13;
            if (f.Mantissa & 0x1000) // Check for rounding
            {
                if (((f.Mantissa & 0x1fff) > 0x1000) || (o.Mantissa & 1)) // if above halfway point or unrounded result is odd
                    o.u++; // Round, might overflow to inf, this is OK
            }
        }
    }

    o.Sign = f.Sign;
    return o;
}

FP16 float_to_half_fast(FP32 f)
{
    FP16 o = { 0 };

    // Based on ISPC reference code (with minor modifications)
    if (f.Exponent == 255) // Inf or NaN (all exponent bits set)
    {
        o.Exponent = 31;
        o.Mantissa = f.Mantissa ? 0x200 : 0; // NaN->qNaN and Inf->Inf
    }
    else // Normalized number
    {
        // Exponent unbias the single, then bias the halfp
        int newexp = f.Exponent - 127 + 15;
        if (newexp >= 31) // Overflow, return signed infinity
            o.Exponent = 31;
        else if (newexp <= 0) // Underflow
        {
            if ((14 - newexp) <= 24) // Mantissa might be non-zero
            {
                uint mant = f.Mantissa | 0x800000; // Hidden 1 bit
                o.Mantissa = mant >> (14 - newexp);
                if ((mant >> (13 - newexp)) & 1) // Check for rounding
                    o.u++; // Round, might overflow into exp bit, but this is OK
            }
        }
        else
        {
            o.Exponent = newexp;
            o.Mantissa = f.Mantissa >> 13;
            if (f.Mantissa & 0x1000) // Check for rounding
                o.u++; // Round, might overflow to inf, this is OK
        }
    }

    o.Sign = f.Sign;
    return o;
}

FP16 float_to_half_fast2(FP32 f)
{
    FP32 infty = { 31 << 23 };
    FP32 magic = { 15 << 23 };
    FP16 o = { 0 };

    uint sign = f.Sign;
    f.Sign = 0;

    // Based on ISPC reference code (with minor modifications)
    if (f.Exponent == 255) // Inf or NaN (all exponent bits set)
    {
        o.Exponent = 31;
        o.Mantissa = f.Mantissa ? 0x200 : 0; // NaN->qNaN and Inf->Inf
    }
    else // (De)normalized number or zero
    {
        f.u &= ~0xfff; // Make sure we don't get sticky bits
        // Not necessarily the best move in terms of accuracy, but matches behavior
        // of other versions.

        // Shift exponent down, denormalize if necessary.
        // NOTE This represents half-float denormals using single precision denormals.
        // The main reason to do this is that there's no shift with per-lane variable
        // shifts in SSE*, which we'd otherwise need. It has some funky side effects
        // though:
        // - This conversion will actually respect the FTZ (Flush To Zero) flag in
        //   MXCSR - if it's set, no half-float denormals will be generated. I'm
        //   honestly not sure whether this is good or bad. It's definitely interesting.
        // - If the underlying HW doesn't support denormals (not an issue with Intel
        //   CPUs, but might be a problem on GPUs or PS3 SPUs), you will always get
        //   flush-to-zero behavior. This is bad, unless you're on a CPU where you don't
        //   care.
        // - Denormals tend to be slow. FP32 denormals are rare in practice outside of things
        //   like recursive filters in DSP - not a typical half-float application. Whether
        //   FP16 denormals are rare in practice, I don't know. Whatever slow path your HW
        //   may or may not have for denormals, this may well hit it.
        f.f *= magic.f;

        f.u += 0x1000; // Rounding bias
        if (f.u > infty.u) f.u = infty.u; // Clamp to signed infinity if overflowed

        o.u = f.u >> 13; // Take the bits!
    }

    o.Sign = sign;
    return o;
}

FP16 float_to_half_fast3(FP32 f)
{
    FP32 f32infty = { 255 << 23 };
    FP32 f16infty = { 31 << 23 };
    FP32 magic = { 15 << 23 };
    uint sign_mask = 0x80000000u;
    uint round_mask = ~0xfffu; 
    FP16 o = { 0 };

    uint sign = f.u & sign_mask;
    f.u ^= sign;

    // NOTE all the integer compares in this function can be safely
    // compiled into signed compares since all operands are below
    // 0x80000000. Important if you want fast straight SSE2 code
    // (since there's no unsigned PCMPGTD).

    if (f.u >= f32infty.u) // Inf or NaN (all exponent bits set)
        o.u = (f.u > f32infty.u) ? 0x7e00 : 0x7c00; // NaN->qNaN and Inf->Inf
    else // (De)normalized number or zero
    {
        f.u &= round_mask;
        f.f *= magic.f;
        f.u -= round_mask;
        if (f.u > f16infty.u) f.u = f16infty.u; // Clamp to signed infinity if overflowed

        o.u = f.u >> 13; // Take the bits!
    }

    o.u |= sign >> 16;
    return o;
}

// Same, but rounding ties to nearest even instead of towards +inf
FP16 float_to_half_fast3_rtne(FP32 f)
{
    FP32 f32infty = { 255 << 23 };
    FP32 f16max   = { (127 + 16) << 23 };
    FP32 denorm_magic = { ((127 - 15) + (23 - 10) + 1) << 23 };
    uint sign_mask = 0x80000000u;
    FP16 o = { 0 };

    uint sign = f.u & sign_mask;
    f.u ^= sign;

    // NOTE all the integer compares in this function can be safely
    // compiled into signed compares since all operands are below
    // 0x80000000. Important if you want fast straight SSE2 code
    // (since there's no unsigned PCMPGTD).

    if (f.u >= f16max.u) // result is Inf or NaN (all exponent bits set)
        o.u = (f.u > f32infty.u) ? 0x7e00 : 0x7c00; // NaN->qNaN and Inf->Inf
    else // (De)normalized number or zero
    {
        if (f.u < (113 << 23)) // resulting FP16 is subnormal or zero
        {
            // use a magic value to align our 10 mantissa bits at the bottom of
            // the float. as long as FP addition is round-to-nearest-even this
            // just works.
            f.f += denorm_magic.f;

            // and one integer subtract of the bias later, we have our final float!
            o.u = f.u - denorm_magic.u;
        }
        else
        {
            uint mant_odd = (f.u >> 13) & 1; // resulting mantissa is odd

            // update exponent, rounding bias part 1
            f.u += ((15 - 127) << 23) + 0xfff;
            // rounding bias part 2
            f.u += mant_odd;
            // take the bits!
            o.u = f.u >> 13;
        }
    }

    o.u |= sign >> 16;
    return o;
}

// Approximate solution. This is faster but converts some sNaNs to
// infinity and doesn't round correctly. Handle with care.
FP16 approx_float_to_half(FP32 f)
{
    FP32 f32infty = { 255 << 23 };
    FP32 f16max = { (127 + 16) << 23 };
    FP32 magic = { 15 << 23 };
    FP32 expinf = { (255 ^ 31) << 23 };
    uint sign_mask = 0x80000000u;
    FP16 o = { 0 };

    uint sign = f.u & sign_mask;
    f.u ^= sign;

    if (!(f.f < f32infty.u)) // Inf or NaN
        o.u = f.u ^ expinf.u;
    else
    {
        if (f.f > f16max.f) f.f = f16max.f;
        f.f *= magic.f;
    }

    o.u = f.u >> 13; // Take the mantissa bits
    o.u |= sign >> 16;
    return o;
}

#ifndef NEON_OPTS
// round-half-up (same as ISPC)
__m128i float_to_half_SSE2(__m128 f)
{
#define CONSTF(name) _mm_castsi128_ps(name)

    __m128i mask_sign       = _mm_set1_epi32(0x80000000u);
    __m128i mask_round      = _mm_set1_epi32(~0xfffu);
    __m128i c_f32infty      = _mm_set1_epi32(255 << 23);
    __m128i c_magic         = _mm_set1_epi32(15 << 23);
    __m128i c_nanbit        = _mm_set1_epi32(0x200);
    __m128i c_infty_as_fp16 = _mm_set1_epi32(0x7c00);
    __m128i c_clamp         = _mm_set1_epi32((31 << 23) - 0x1000);

    __m128  msign       = CONSTF(mask_sign);
    __m128  justsign    = _mm_and_ps(msign, f);
    __m128i f32infty    = c_f32infty;
    __m128  absf        = _mm_xor_ps(f, justsign);
    __m128  mround      = CONSTF(mask_round);
    __m128i absf_int    = _mm_castps_si128(absf); // pseudo-op, but val needs to be copied once so count as mov
    __m128i b_isnan     = _mm_cmpgt_epi32(absf_int, f32infty);
    __m128i b_isnormal  = _mm_cmpgt_epi32(f32infty, _mm_castps_si128(absf));
    __m128i nanbit      = _mm_and_si128(b_isnan, c_nanbit);
    __m128i inf_or_nan  = _mm_or_si128(nanbit, c_infty_as_fp16);

    __m128  fnosticky   = _mm_and_ps(absf, mround);
    __m128  scaled      = _mm_mul_ps(fnosticky, CONSTF(c_magic));
    __m128  clamped     = _mm_min_ps(scaled, CONSTF(c_clamp)); // logically, we want PMINSD on "biased", but this should gen better code
    __m128i biased      = _mm_sub_epi32(_mm_castps_si128(clamped), _mm_castps_si128(mround));
    __m128i shifted     = _mm_srli_epi32(biased, 13);
    __m128i normal      = _mm_and_si128(shifted, b_isnormal);
    __m128i not_normal  = _mm_andnot_si128(b_isnormal, inf_or_nan);
    __m128i joined      = _mm_or_si128(normal, not_normal);

    __m128i sign_shift  = _mm_srli_epi32(_mm_castps_si128(justsign), 16);
    __m128i final       = _mm_or_si128(joined, sign_shift);

    // ~20 SSE2 ops
    return final;

#undef CONSTF
}

// round-to-nearest-even
// this is an adaptation of float_to_half_fast3_rtne which is the code
// you should read to understand the algorithm.
__m128i float_to_half_rtne_SSE2(__m128 f)
{
#define CONSTF(name) _mm_castsi128_ps(name)

    __m128i mask_sign       = _mm_set1_epi32(0x80000000u);
    __m128i c_f16max        = _mm_set1_epi32((127 + 16) << 23); // all FP32 values >=this round to +inf
    __m128i c_nanbit        = _mm_set1_epi32(0x200);
    __m128i c_infty_as_fp16 = _mm_set1_epi32(0x7c00);
    __m128i c_min_normal    = _mm_set1_epi32((127 - 14) << 23); // smallest FP32 that yields a normalized FP16
    __m128i c_subnorm_magic = _mm_set1_epi32(((127 - 15) + (23 - 10) + 1) << 23);
    __m128i c_normal_bias   = _mm_set1_epi32(0xfff - ((127 - 15) << 23)); // adjust exponent and add mantissa rounding

    __m128  msign       = CONSTF(mask_sign);
    __m128  justsign    = _mm_and_ps(msign, f);
    __m128  absf        = _mm_xor_ps(f, justsign);
    __m128i absf_int    = _mm_castps_si128(absf); // the cast is "free" (extra bypass latency, but no thruput hit)
    __m128i f16max      = c_f16max;
    __m128  b_isnan     = _mm_cmpunord_ps(absf, absf); // is this a NaN?
    __m128i b_isregular = _mm_cmpgt_epi32(f16max, absf_int); // (sub)normalized or special?
    __m128i nanbit      = _mm_and_si128(_mm_castps_si128(b_isnan), c_nanbit);
    __m128i inf_or_nan  = _mm_or_si128(nanbit, c_infty_as_fp16); // output for specials

    __m128i min_normal  = c_min_normal;
    __m128i b_issub     = _mm_cmpgt_epi32(min_normal, absf_int);

    // "result is subnormal" path
    __m128  subnorm1    = _mm_add_ps(absf, CONSTF(c_subnorm_magic)); // magic value to round output mantissa
    __m128i subnorm2    = _mm_sub_epi32(_mm_castps_si128(subnorm1), c_subnorm_magic); // subtract out bias

    // "result is normal" path
    __m128i mantoddbit  = _mm_slli_epi32(absf_int, 31 - 13); // shift bit 13 (mantissa LSB) to sign
    __m128i mantodd     = _mm_srai_epi32(mantoddbit, 31); // -1 if FP16 mantissa odd, else 0

    __m128i round1      = _mm_add_epi32(absf_int, c_normal_bias);
    __m128i round2      = _mm_sub_epi32(round1, mantodd); // if mantissa LSB odd, bias towards rounding up (RTNE)
    __m128i normal      = _mm_srli_epi32(round2, 13); // rounded result

    // combine the two non-specials
    __m128i nonspecial  = _mm_or_si128(_mm_and_si128(subnorm2, b_issub), _mm_andnot_si128(b_issub, normal));

    // merge in specials as well
    __m128i joined      = _mm_or_si128(_mm_and_si128(nonspecial, b_isregular), _mm_andnot_si128(b_isregular, inf_or_nan));

    __m128i sign_shift  = _mm_srli_epi32(_mm_castps_si128(justsign), 16);
    __m128i final       = _mm_or_si128(joined, sign_shift);

    // ~28 SSE2 ops
    return final;

#undef CONSTF
}

__m128i approx_float_to_half_SSE2(__m128 f)
{
#define CONSTF(name) _mm_castsi128_ps(name)

    __m128i mask_fabs  = _mm_set1_epi32(0x7fffffff);
    __m128i c_f32infty = _mm_set1_epi32((255 << 23));
    __m128i c_expinf   = _mm_set1_epi32((255 ^ 31) << 23);
    __m128i c_f16max   = _mm_set1_epi32((127 + 16) << 23);
    __m128i c_magic    = _mm_set1_epi32(15 << 23);

    __m128  mabs        = CONSTF(mask_fabs);
    __m128  fabs        = _mm_and_ps(mabs, f);
    __m128  justsign    = _mm_xor_ps(f, fabs);

    __m128  f16max      = CONSTF(c_f16max);
    __m128  expinf      = CONSTF(c_expinf);
    __m128  infnancase  = _mm_xor_ps(expinf, fabs);
    __m128  clamped     = _mm_min_ps(f16max, fabs);
    __m128  b_notnormal = _mm_cmpnlt_ps(fabs, CONSTF(c_f32infty));
    __m128  scaled      = _mm_mul_ps(clamped, CONSTF(c_magic));

    __m128  merge1      = _mm_and_ps(infnancase, b_notnormal);
    __m128  merge2      = _mm_andnot_ps(b_notnormal, scaled);
    __m128  merged      = _mm_or_ps(merge1, merge2);

    __m128i shifted     = _mm_srli_epi32(_mm_castps_si128(merged), 13);
    __m128i signshifted = _mm_srli_epi32(_mm_castps_si128(justsign), 16);
    __m128i final       = _mm_or_si128(shifted, signshifted);

    // ~15 SSE2 ops
    return final;

#undef CONSTF
}
#endif

// from fox toolkit float->half code (which "approx" variants match)
static uint basetable[512];
static unsigned char shifttable[512];

void fp16_generatetables()
{
    unsigned int i;
    int e;
    for(i=0; i<256; ++i){
        e=i-127;
        if(e<-24){                  // Very small numbers map to zero
            basetable[i|0x000]=0x0000;
            basetable[i|0x100]=0x8000;
            shifttable[i|0x000]=24;
            shifttable[i|0x100]=24;
        }
        else if(e<-14){             // Small numbers map to denorms
            basetable[i|0x000]=(0x0400>>(-e-14));
            basetable[i|0x100]=(0x0400>>(-e-14)) | 0x8000;
            shifttable[i|0x000]=-e-1;
            shifttable[i|0x100]=-e-1;
        }
        else if(e<=15){             // Normal numbers just lose precision
            basetable[i|0x000]=((e+15)<<10);
            basetable[i|0x100]=((e+15)<<10) | 0x8000;
            shifttable[i|0x000]=13;
            shifttable[i|0x100]=13;
        }
        else if(e<128){             // Large numbers map to Infinity
            basetable[i|0x000]=0x7C00;
            basetable[i|0x100]=0xFC00;
            shifttable[i|0x000]=24;
            shifttable[i|0x100]=24;
        }
        else{                       // Infinity and NaN's stay Infinity and NaN's
            basetable[i|0x000]=0x7C00;
            basetable[i|0x100]=0xFC00;
            shifttable[i|0x000]=13;
            shifttable[i|0x100]=13;
        }
    }
}

// also from fox toolkit
uint float_to_half_foxtk(uint f)
{
    return basetable[(f>>23)&0x1ff]+((f&0x007fffff)>>shifttable[(f>>23)&0x1ff]);
}

// from half->float code - just for verification.
FP32 half_to_float(FP16 h)
{
    static const FP32 magic = { 113 << 23 };
    static const uint shifted_exp = 0x7c00 << 13; // exponent mask after shift
    FP32 o;

    o.u = (h.u & 0x7fff) << 13;     // exponent/mantissa bits
    uint exp = shifted_exp & o.u;   // just the exponent
    o.u += (127 - 15) << 23;        // exponent adjust

    // handle exponent special cases
    if (exp == shifted_exp) // Inf/NaN?
        o.u += (128 - 16) << 23;    // extra exp adjust
    else if (exp == 0) // Zero/Denormal?
    {
        o.u += 1 << 23;             // extra exp adjust
        o.f -= magic.f;             // renormalize
    }

    o.u |= (h.u & 0x8000) << 16;    // sign bit
    return o;
}

FP32 half_to_float_lit(unsigned short u)
{
    FP16 fp16 = { u };
    return half_to_float(fp16);
}

#ifdef FP16_MAIN
int main(int argc, char **argv)
{
    FP32 f;
    FP16 full, fast, fast2, fast3;
    uint u;

    generatetables();

#if 0 // commented out since one full pass is slow enough...
    u = 0;

    //u = 0x32000000;
    //u = 0x32fff800;
    //u = 0x33000000;

    do
    {
        f.u = u;
        full = float_to_half_full(f);
        fast = float_to_half_fast(f);
        fast2 = float_to_half_fast2(f);
        fast3 = float_to_half_fast3(f);
        if (full.u != fast.u || full.u != fast2.u || full.u != fast3.u)
        {
            FP32 fastback, fast2back, fast3back;
            fastback = half_to_float(fast);
            fast2back = half_to_float(fast2);
            fast3back = half_to_float(fast3);
            printf("mismatch! val=%08x full=%04x fast=%04x->%08x fast2=%04x->%08x fast3=%04x->%08x\n", u, full.u,
                fast.u, fastback.u, fast2.u, fast2back.u, fast3.u, fast3back.u);
            return 1;
        }

        ++u;
        if ((u & 0xffffff) == 0)
            printf("  %02x\n", (u-1) >> 24);
    }
    while (u);
#endif

#if 0
    printf("ISPC vs. round-to-nearest-even:\n");
    int num_expected_mismatch = 0;
    u = 0;
    do
    {
        f.u = u;
        full = float_to_half_full(f);
        FP16 rtne = float_to_half_full_rtne(f);

        if (full.u != rtne.u)
        {
            // Some mismatches expected, but make sure this is one of them!
            FP32 f_full = half_to_float(full);
            FP32 f_rtne = half_to_float(rtne);

            // Expected cases: f_full and f_rtne are equidistant from true value (in ulps), full is odd and rtne is even.
            // (Unless we rounded to +-0.)
            int diff_full = abs((int) (f_full.u - u));
            int diff_rtne = abs((int) (f_rtne.u - u));
            if ((diff_full == diff_rtne || (diff_full == 0x800000 && (rtne.u & 0x7fff) == 0)) && ((full.u & 1) == 1) && ((rtne.u & 1) == 0))
                ++num_expected_mismatch;
            else
            {
                printf("unexpected mismatch! val=%08x mant=%06x full=%04x rtne=%04x diff_full=%x diff_rtne=%x\n", u, f.Mantissa, full.u, rtne.u, diff_full, diff_rtne);
                return 1;
            }
        }

        ++u;
        if ((u & 0xffffff) == 0)
            printf("  %02x\n", (u-1) >> 24);
    } while (u);
    printf("%d expected mismatches\n", num_expected_mismatch);
#endif

#if 0 // RTNE vs. ref
    u = 0;

    do
    {
        f.u = u;
        full = float_to_half_full_rtne(f);
        fast3 = float_to_half_fast3_rtne(f);
        if (full.u != fast3.u)
        {
            FP32 fast3back = half_to_float(fast3);
            printf("mismatch! val=%08x full=%04x fast3=%04x->%08x\n", u, full.u, fast3.u, fast3back.u);
            return 1;
        }

        ++u;
        if ((u & 0xffffff) == 0)
            printf("  %02x\n", (u-1) >> 24);
    }
    while (u);
#endif

#if 0
    printf("SSE2:\n");
    u = 0;
    do
    {
        __m128 ssein;
        __m128i ref, sseout;

        for (int j=0; j < 4; j++)
        {
            ssein.m128_u32[j] = u + j;
            f.u = u + j;
            full = float_to_half_full(f);
            ref.m128i_u32[j] = full.u;
        }

        sseout = float_to_half_SSE2(ssein);

        for (int j=0; j < 4; j++)
        {
            if (sseout.m128i_u32[j] != ref.m128i_u32[j])
            {
                printf("mismatch! val=%08x full=%04x fastSSE2=%04x\n", u+j, ref.m128i_u32[j], sseout.m128i_u32[j]);
                return 1;
            }
        }

        u += 4;
        if ((u & 0xffffff) == 0)
            printf("  %02x\n", (u-1) >> 24);
    } while (u);
#endif

#if 0
    printf("approx:\n");
    u = 0;
    do
    {
        __m128 ssein;
        __m128i ref, sseout;

        for (int j=0; j < 4; j++)
        {
            uint x = u + j;
            ssein.m128_u32[j] = x;
            ref.m128i_u32[j] = float_to_half_foxtk(x);
        }

        sseout = approx_float_to_half_SSE2(ssein);

        for (int j=0; j < 4; j++)
        {
            if (abs((int) (sseout.m128i_u32[j] - ref.m128i_u32[j])) > 1)
            {
                printf("mismatch! val=%08x ref=%04x approx=%04x\n", ssein.m128_u32[j], ref.m128i_u32[j], sseout.m128i_u32[j]);
                printf("exp = %d\n", ((ssein.m128_u32[j] >> 23) & 255) - 127);
                return 1;
            }
        }

        u += 4;
        if ((u & 0xffffff) == 0)
            printf("  %02x\n", (u - 1) >> 24);

    } while (u);
#endif

#if 1
    printf("SSE2 RTNE:\n");
    u = 0;
    do
    {
        __m128 ssein;
        __m128i ref, sseout;

        for (int j=0; j < 4; j++)
        {
            ssein.m128_u32[j] = u + j;
            f.u = u + j;
            full = float_to_half_full_rtne(f);
            ref.m128i_u32[j] = full.u;
        }

        sseout = float_to_half_rtne_SSE2(ssein);

        for (int j=0; j < 4; j++)
        {
            if (sseout.m128i_u32[j] != ref.m128i_u32[j])
            {
                printf("mismatch! val=%08x full=%04x fastSSE2=%04x\n", u+j, ref.m128i_u32[j], sseout.m128i_u32[j]);
                return 1;
            }
        }

        u += 4;
        if ((u & 0xffffff) == 0)
            printf("  %02x\n", (u-1) >> 24);
    } while (u);
#endif

#if 0
    // The "Steinar test" :)
    {
        FP32 val1 = half_to_float_lit(0x3c00);
        FP32 val2 = half_to_float_lit(0x3c01);

        FP32 sum;
        sum.f = val1.f + val2.f;

        printf("sum rtne:  0x%04x (should be 0x4000)\n", float_to_half_full_rtne(sum).u);
        
        FP32 tiny;
        tiny.f = 0.5f*5.9604644775390625e-08f;
        printf("tiny rtne: 0x%04x (should be 0x0000)\n", float_to_half_full_rtne(tiny).u);
    }
#endif

    return 0;
}
#endif
