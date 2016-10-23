/*
This file is part of libcsdr.

	Copyright (c) Andras Retzler, HA7ILM <randras@sdr.hu>
	Copyright (c) Warren Pratt, NR0V <warren@wpratt.com>
	Copyright 2006,2010,2012 Free Software Foundation, Inc.

    libcsdr is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    libcsdr is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with libcsdr.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#ifdef LIBCSDR_GPL

typedef struct shift_addition_data_s
{
	float sindelta;
	float cosdelta;
	float rate;
} shift_addition_data_t;
shift_addition_data_t shift_addition_init(float rate);
float shift_addition_cc(complexf *input, complexf* output, int input_size, shift_addition_data_t d, float starting_phase);
float shift_addition_fc(float *input, complexf* output, int input_size, shift_addition_data_t d, float starting_phase);
void shift_addition_cc_test(shift_addition_data_t d);

float agc_ff(float* input, float* output, int input_size, float reference, float attack_rate, float decay_rate, float max_gain, short hang_time, short attack_wait_time, float gain_filter_alpha, float last_gain);

typedef struct decimating_shift_addition_status_s
{
	int decimation_remain;
	float starting_phase;
	int output_size;
} decimating_shift_addition_status_t;
decimating_shift_addition_status_t decimating_shift_addition_cc(complexf *input, complexf* output, int input_size, shift_addition_data_t d, int decimation, decimating_shift_addition_status_t s);
shift_addition_data_t decimating_shift_addition_init(float rate, int decimation);

#endif

