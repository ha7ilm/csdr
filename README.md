libcsdr
=======

*libcsdr* is a set of simple DSP routines for Software Defined Radio.  
It is mostly useful for AM/FM/SSB demodulation and spectrum display.  
Feel free to use it in your projects.  
Most of the code is available under the permissive BSD license, with some optional parts under GPL. For additional details, see <a href="#licensing">licensing</a>.

- The package comes with a command-line tool `csdr`, which lets you build DSP processing chains by shell pipes.
- The code of *libcsdr* was intended to be easy to follow.
- *libcsdr* was designed to use auto-vectorization available in *gcc*. It means that it can achieve some speedup by taking advantage of SIMD command sets available in today's CPUs (e.g. SSE on x86 and NEON on ARM).

How to compile
--------------
The project was only tested on Linux. It has the following dependencies: `libfftw3-dev`

	make
	sudo make install

If you compile on ARM, please edit the Makefile and tailor `PARAMS_NEON` for your CPU.
To run the examples, you will also need <a href="http://sdr.osmocom.org/trac/wiki/rtl-sdr">rtl_sdr</a> from Osmocom, and the following packages (at least on Debian): `mplayer octave gnuplot gnuplot-x11`

Credits
-------
The library was written by Andras Retzler, HA7ILM &lt;<randras@sdr.hu>&gt;.

I would like to say special thanks to Péter Horváth, HA5CQA and János Selmeczi, HA5FT for their continous help and support.

Usage by example
----------------

### Demodulate WFM

	rtl_sdr -s 240000 -f 89500000 -g 20 - | csdr convert_u8_f | csdr fmdemod_quadri_cf | csdr fractional_decimator_ff 5 | csdr deemphasis_wfm_ff 48000 50e-6 | csdr convert_f_i16 | mplayer -cache 1024 -quiet -rawaudio samplesize=2:channels=1:rate=48000 -demuxer rawaudio -

- Baseband I/Q signal is coming from an RTL-SDR USB dongle, with a center frequency of `-f 104300000` Hz, a sampling rate of `-s 240000` samples per second.
- The `rtl_sdr` tool outputs an unsigned 8-bit I/Q signal (one byte of I sample and one byte of Q coming after each other), but `libcsdr` DSP routines internally use floating point data type, so we convert the data stream of `unsigned char` to `float` by `csdr convert_u8_f`.
- We want to listen one radio station at the frequency `-f 89500000` Hz (89.5 MHz).
- No other radio station is within the sampled bandwidth, so we send the signal directly to the demodulator. (This is an easy, but not perfect solution as the anti-aliasing filter at RTL-SDR DDC is too short.) 
- After FM demodulation we decimate the signal by a factor of 5 to match the rate of the audio card (240000 / 5 = 48000).
- A de-emphasis filter is used, because pre-emphasis is applied at the transmitter to compensate noise at higher frequencies. The time constant for de-emphasis for FM broadcasting in Europe is 50 microseconds (hence the `50e-6`).
- Also, `mplayer` cannot play floating point audio, so we convert our signal to a stream of 16-bit integers.  

### Demodulate WFM: advanced

	rtl_sdr -s 2400000 -f 89300000 -g 20 - | csdr convert_u8_f | csdr shift_addition_cc -0.17 | csdr fir_decimate_cc 10 0.05 HAMMING | csdr fmdemod_quadri_cf | csdr fractional_decimator_ff 5 | csdr deemphasis_wfm_ff 48000 50e-6 | csdr convert_f_i16 | mplayer -cache 1024 -quiet -rawaudio samplesize=2:channels=1:rate=48000 -demuxer rawaudio -

- We want to listen to one radio station, but input signal contains multiple stations, and its bandwidth is too large for sending it directly to the FM demodulator.
- We shift the signal to the center frequency of the station we want to receive: `-0.17*2400000*0.5 = -204000`, so basically we will listen to the radio station centered at 89504000 Hz.
- We decimate the signal by a factor of 10. The rolloff for the FIR filter used for decimation will be 10% of total bandwidth (as of parameter 0.05 is 10% of 0.5). Hamming window will be used for windowed FIR filter design.

Sample rates look like this:


	             2.4 Msps                     240 ksps                                  48 ksps
	I/Q source ------------> FIR decimation ------------> FM demod -> frac. decimation ---------> deemphasis -> sound card


*Note:* there is an example shell script that does this for you (without the unnecessary shift operation). If you just want to listen to FM radio, type:

	csdr-fm 89.5 20 

The first parameter is the frequency in MHz, and the second optional parameter is the RTL-SDR tuner gain in dB.

### Demodulate NFM
	
	rtl_sdr -s 2400000 -f 145000000 -g 20 - | csdr convert_u8_f | csdr shift_addition_cc `python -c "print (145000000-145350000)/(0.5*2400000)"` | csdr fir_decimate_cc 50 0.005 HAMMING | csdr fmdemod_quadri_cf | csdr limit_ff | csdr deemphasis_nfm_ff 48000 | csdr fastagc_ff | csdr convert_f_i16 | mplayer -cache 1024 -quiet -rawaudio samplesize=2:channels=1:rate=48000 -demuxer rawaudio -

- Note that the decimation factor is higher (we want to select a ~25 kHz channel).
- Also there is a python hack to calculate the relative shift offset. The real receiver frequency is `145350000` Hz.
- The de-emphasis filter is a fixed FIR filter that has a passband of 400-4000 Hz, also with a roll-off of -20 dB/decade.

### Demodulate AM

	rtl_sdr -s 2400000 -f 145000000 -g 20 - | csdr convert_u8_f | csdr shift_addition_cc `python -c "print (145000000-144400000)/(0.5*2400000)"` | csdr fir_decimate_cc 50 0.005 HAMMING | csdr amdemod_cf | csdr fastdcblock_ff | csdr agc_ff | csdr limit_ff | csdr convert_f_i16 | mplayer -cache 1024 -quiet -rawaudio samplesize=2:channels=1:rate=48000 -demuxer rawaudio -

- `amdemod_cf` is used as demodulator. 
- `agc_ff` should be used for AM and SSB.

### Design FIR band-pass filter (with complex taps)

	csdr firdes_bandpass_c 0 0.5 59 HAMMING --octave | octave -i

- ...and then plot its frequency response with octave. 
- It will design a filter that lets only the positive frequencies pass (low cut is 0, high cut is 0.5 - these are relative to the sampling rate).
- If `--octave` and everything that follows is removed from the command, you get only the taps. E. g. the raw output of `firdes_lowpass_f` can be easily copied to C code.

### Demodulate SSB

	rtl_sdr -s 2400000 -f 145000000 -g 20 - | csdr convert_u8_f | csdr shift_addition_cc `python -c "print (145000000-144400000)/(0.5*2400000)"` | csdr fir_decimate_cc 50 0.005 HAMMING | csdr bandpass_fir_fft_cc 0 0.1 0.05 | csdr realpart_cf | csdr agc_ff | csdr limit_ff | csdr convert_f_i16 | mplayer -cache 1024 -quiet -rawaudio samplesize=2:channels=1:rate=48000 -demuxer rawaudio -

- It is a modified Weaver-demodulator. The complex FIR filter removes the lower sideband and lets only the upper pass (USB). If you want to demodulate LSB, change `bandpass_fir_fft_cc 0 0.05` to `bandpass_fir_fft_cc -0.05 0`. 

### Draw FFT

	rtl_sdr -s 2400000 -f 104300000 -g 20 - | csdr convert_u8_f | csdr fft_cc 1024 1200000 HAMMING --octave | octave -i > /dev/null

- We calculate the Fast Fourier Transform by `csdr fft_cc` on the first 1024 samples of every block of 1200000 complex samples coming after each other. (We calculate FFT from 1024 samples and then skip 1200000-1024=1198976 samples. This way we will calculate FFT two times every second.)
- The window used for FFT is the Hamming window, and the output consists of commands that can be directly interpreted by GNU Octave which plots us the spectrum. 

Usage
-----
Some basic concepts on using *libcsdr*:

### Data types
Function name endings found in *libcsdr* mean the input and output data types of the particular function. (This is similar to GNU Radio naming conventions).
Data types are noted as it follows:

- `f` is `float` (single percision)
- `c` is `complexf` (two single precision floating point values in a struct) 
- `u8` is `unsigned char` of 1 byte/8 bits (e. g. the output of `rtl_sdr` is of `u8`)
- `i16` is `signed short` of 2 bytes/16 bits (e. g. sound card input is usually `i16`)

Functions usually end as:

- `_ff` float input, float output
- `_cf` complex input, float output
- `_cc` complex input, complex output

Regarding *csdr*, it can convert a real/complex stream from one data format to another, to interface it with other SDR tools and the sound card. 
The following commands are available: 

- `csdr convert_u8_f` 
- `csdr convert_f_u8` 
- `csdr convert_i16_f` 
- `csdr convert_f_i16`

How to interpret: `csdr convert_<src>_<dst>`
You can use these commands on complex streams, too, as they are only interleaved values (I,Q,I,Q,I,Q... coming after each other).

#### csdr commands

`csdr` should be considered as a reference implementation on using `libcsdr`. For additional details on how to use the library, check `csdr.c` and `libcsdr.c`.

Regarding `csdr`, the first command-line parameter is the name of a function, others are the parameters for the given function. Compulsory parameters are noted as `<parameter>`, optional parameters are noted as `[parameter]`.
Optional parameters have safe defaults, for more info look at the code.

	realpart_cf

It takes the real part of the complex signal, and throws away the imaginary part.

	clipdetect_ff

It clones the signal (the input and the output is the same), but it prints a warning on `stderr` if any sample value is out of the -1.0 ... 1.0 range.

	limit_ff [max_amplitude]

The input signal amplitude will not be let out of the `-max_amplitude ... max_amplitude` range.

	gain_ff <gain>

It multiplies all samples by `gain`.

	clone

It copies the input to the output.

	yes_f <to_repeat> [buf_times]

It outputs continously the `to_repeat` float number. 
If `buf_times` is not given, it never stops.
Else, after outputing `buf_times` number of buffers (the size of which is stated in the `BUFSIZE` macro), it exits.

	shift_math_cc <rate>

It shifts the complex spectrum by `rate`.
`rate` is a floating point number between -0.5 and 0.5. 
`rate` is relative to the sampling rate. 

Internally, a sine and cosine wave is generated to perform this function, and this function uses `math.h` for this purpose, which is quite accurate, but not always very fast.

	shift_addition_cc <rate>

Operation is the same as with `shift_math_cc`.

Internally, this function uses trigonometric addition formulas to generate sine and cosine, which is a bit faster. (About 4 times on the machine I have tested it on.)

	shift_addition_cc_test

This function was used to test the accuracy of the method above.

	dcblock_ff

This is a DC blocking IIR filter.

	fastdcblock_ff

This is a DC blocker that works based on the average of the buffer.

	fmdemod_atan_cf

It is an FM demodulator that internally uses the `atan` function in `math.h`, so it is not so fast.

	fmdemod_quadri_cf

It is an FM demodulator that is based on the quadri-correlator method, and it can be effectively auto-vectorized, so it should be faster.

	fmdemod_quadri_novect_cf

It has more easily understandable code than the previous one, but can't be auto-vectorized.

	deemphasis_wfm_ff <sample_rate> <tau>

It does de-emphasis with the given RC time constant `tau`.
Different parts of the world use different pre-emphasis filters for FM broadcasting.
In Europe, `tau` should be chosen as `50e-6`, and in the USA, `tau` should be `75e-6`.

	deemphasis_nfm_ff <one_of_the_predefined_sample_rates>

It does de-emphasis on narrow-band FM for communication equipment (e.g. two-way radios).
It uses fixed filters so it works only on predefined sample rates, for the actual list of them run: `cat libcsdr.c | grep DNFMFF_ADD_ARRAY`

	amdemod_cf

It is an AM demodulator that uses `sqrt`. On some architectures `sqrt` can be directly calculated by dedicated CPU instructions, but on others it may be slower. 

amdemod_estimator_cf

It is an AM demodulator that uses an estimation method that is faster but less accurate than `amdemod_cf`.

	firdes_lowpass_f <cutoff_rate> <length> [window [--octave]]

Low-pass FIR filter design function to output real taps, with a `cutoff_rate` proportional to the sampling frequency, using the windowed sinc filter design method.
`cutoff_rate` can be between 0 and 0.5.

`length` is the number of filter taps to output, and should be odd.
The longer the filter kernel is, the shorter the transition bandwidth is, but the more CPU time it takes to process the filter.
The transition bandwidth (proportional to the sampling rate) can be calculated as: `transition_bw = 4 / length`.
Some functions (below) require the `transition_bw` to be given instead of filter `length`. Try to find the best compromise between speed and accuracy by changing this parameter.

`window` is the window function used to compensate finite filter length. Its typical values are: `HAMMING`, `BLACKMAN`, `BOXCAR`. For the actual list of values, run: `cpp libcsdr.c | grep window\ ==`

The `--octave` parameter lets you directly view the filter response in `octave`. For more information, look at the [Usage by example] section.

	firdes_bandpass_c <low_cut> <high_cut> <length> [window [--octave]]

Band-pass FIR filter design function to output complex taps.
`low_cut` and ` high_cut` both may be between -0.5 and 0.5, and are also proportional to the sampling frequency.

Other parameters were explained above at `firdes_lowpass_f`.

	fir_decimate_cc <decimation_factor> [transition_bw [window]]

It is a decimator that keeps one sample out of `decimation_factor` samples. 
To avoid aliasing, it runs a filter on the signal and removes spectral components above `0.5 × nyquist_frequency × decimation_factor`.

`transition_bw` and `window` are the parameters of the filter.

	rational_resampler_ff <interpolation> <decimation> [transition_bw [window]]

It is a resampler that takes integer values of `interpolation` and `decimation`.
The output sample rate will be `interpolation / decimation × input_sample_rate`.

`transition_bw` and `window` are the parameters of the filter.

	fractional_decimator_ff <decimation_rate> [transition_bw [window]]

It can decimate by a floating point ratio.

`transition_bw` and `window` are the parameters of the filter.

	bandpass_fir_fft_cc <low_cut> <high_cut> <transition_bw> [window]

It performs a bandpass FIR filter on complex samples, using FFT and the overlap-add method.

Parameters are described under `firdes_bandpass_c` and `firdes_lowpass_f`.

	agc_ff [hang_time [reference [attack_rate [decay_rate [max_gain [attack_wait [filter_alpha]]]]]]]

It is an automatic gain control function. 

- `hang_time` is the number of samples to wait before strating to increase the gain after a peak.
- `reference` is the reference level for the AGC. It tries to keep the amplitude of the output signal close to that.
- `attack_rate` is the rate of decreasing the signal level if it gets higher than it used to be before.
- `decay_rate` is the rate of increasing the signal level if it gets lower than it used to be before.
- AGC won't increase the gain over `max_gain`.
- `attack_wait` is the number of sampels to wait before starting to decrease the gain, because sometimes very short peaks happen, and we don't want them to spoil the reception by substantially decreasing the gain of the AGC.
- `filter_alpha` is the parameter of the loop filter.

Its default parameters work best for an audio signal sampled at 48000 Hz.

	fastagc_ff [block_size [reference]]

It is a faster AGC that linearly changes the gain, taking the highest amplitude peak in the buffer into consideration. Its output will never exceed `-reference ... reference`.

	fft_cc <fft_size> <out_of_every_n_samples> [window [--octave] [--benchmark]]

It performs an FFT on the first `fft_size` samples out of `out_of_every_n_samples`, thus skipping `out_of_every_n_samples - fft_size` samples in the input. 

It can draw the spectrum by using `--octave`, for more information, look at the [Usage by example] section.

FFTW can be faster if we let it optimalize a while before starting the first transform, hence the `--benchmark` switch.

	fft_benchmark <fft_size> <fft_cycles> [--benchmark]

It measures the time taken to process `fft_cycles` transforms of `fft_size`.
It lets FFTW optimalize if used with the `--benchmark` switch.

	lowpower_cf [add_db]

Calculates `10*log10(i^2+q^2)+add_db` for the input complex samples. It is useful for drawing power spectrum graphs.

#### Control via pipes

Some parameters can be changed while the `csdr` process is running. To achieve this, some `csdr` functions have special parameters. You have to supply a fifo previously created by the `mkfifo` command. Processing will only start after the first control command has been received by `csdr` over the FIFO.

	shift_addition_cc --fifo <fifo_path>

By writing to the given FIFO file with the syntax below, you can control the shift rate:

	<shift_rate>\n

E.g. you can send `-0.3\n`

Processing will only start after the first control command has been received by `csdr` over the FIFO.

	bandpass_fir_fft_cc --fifo <fifo_path> <transition_bw> [window]

By writing to the given FIFO file with the syntax below, you can control the shift rate:

	<low_cut> <high_cut>\n
	
E.g. you can send `-0.05 0.02\n`

#### Testbench

`csdr` was tested with GNU Radio Companion flowgraphs. These flowgraphs are available under the directory `grc_tests`, and they require the `gr-ha5kfu` set of blocks for GNU Radio.  

## [Licensing] (#licensing)

Most of the code for `libcsdr` is under BSD license.  
However, before the implementation of some algoritms, GPL-licensed code from other applications have been reviewed.
In order to eliminate any licesing issues, these parts are placed under a different file.
However, the library is still fully functional with BSD-only code, altough having only less-optimized versions of some algorithms.  
It should also be noted that if you compile with `-DUSE_FFTW` and `-DLIBCSDR_GPL` (as default), the GPL license would apply on the whole result.

