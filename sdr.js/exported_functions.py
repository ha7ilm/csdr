"""
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
"""


exported_functions= \
"""firdes_lowpass_f
firdes_bandpass_c
firdes_wkernel_blackman
firdes_wkernel_hamming
firdes_wkernel_boxcar
firdes_get_window_from_string
firdes_get_string_from_window
firdes_filter_len
fmdemod_quadri_cf
fmdemod_quadri_novect_cf
fmdemod_atan_cf
amdemod_cf
amdemod_estimator_cf
limit_ff
fir_decimate_cc
deemphasis_nfm_ff
deemphasis_wfm_ff
shift_math_cc
dcblock_ff
fastdcblock_ff
fastagc_ff
rational_resampler_ff
rational_resampler_get_lowpass_f
apply_window_c
apply_window_f
logpower_cf
fractional_decimator_ff
shift_table_deinit
shift_table_init
shift_table_cc
log2n
next_pow2
apply_fir_fft_cc
gain_ff
convert_u8_f
convert_f_u8
convert_f_i16
convert_i16_f
shift_addition_init
shift_addition_cc
shift_addition_cc_test
agc_ff
decimating_shift_addition_cc
decimating_shift_addition_init
encode_ima_adpcm_i16_u8
decode_ima_adpcm_u8_i16"""

exported_functions_quoted=map(lambda x:"'_"+x+"'",exported_functions.split("\n"))
print "["+(", ".join(exported_functions_quoted))+"]"
