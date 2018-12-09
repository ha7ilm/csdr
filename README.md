CSDR
====

`csdr` is a command line tool to carry out DSP tasks for Software Defined Radio.

It can be used to build simple signal processing flow graphs, right from the command line. 

The included `libcsdr` library contains the DSP functions that `csdr` makes use of. It was designed to use auto-vectorization available in `gcc`, and also has some functions optimized with inline assembly for ARM NEON to achieve some speedup by taking advantage of SIMD command sets available in today's CPUs.

Feel free to use it in your projects.  
Most of the code is available under the permissive BSD license, with some optional parts under GPL. For additional details, see <a href="#licensing">licensing</a>.

`csdr` has already been used to build:

- AM, FM, SSB, CW and BPSK31 demodulators and waterfall display in [OpenWebRX](https://github.com/simonyiszk/openwebrx),
- AM, FM, SSB modulators in [qtcsdr](https://github.com/ha7ilm/qtcsdr) that can also be used standalone with [rpitx](https://github.com/ha7ilm/rpitx-app-note),
- a demodulator for FSK transmissions sent with the CC1111 wireless MCU, and also a standalone RTTY demodulator. 

This animation shows the Gardner timing recovery algorithm in `csdr` locking on a baseband BPSK signal:

![Gardner](https://raw.githubusercontent.com/wiki/simonyiszk/csdr/gardner.gif)

(The symbol is sampled at the left red dot. The algorithm moves the middle dot as close to the symbol transition center, as possible.)

How to compile
--------------

    make
    sudo make install

The project was only tested on Linux. It has the following dependencies: `libfftw3-dev`

If you compile on ARM, please edit the Makefile and tailor `PARAMS_NEON` for your CPU.

To run the examples, you will also need <a href="http://sdr.osmocom.org/trac/wiki/rtl-sdr">rtl_sdr</a> from Osmocom, and the following packages (at least on Debian): `mplayer octave gnuplot gnuplot-x11`

If you compile `fftw3` from sources for use with `libcsdr`, you need to configure it with 32-bit float support enabled: 

    ./configure --enable-float

(This is for `fftw3`, not `libcsdr`. You do not need to run the configure script before compiling `libcsdr`.)

Credits
-------
The library was written by Andras Retzler, HA7ILM &lt;<randras@sdr.hu>&gt;.

I would like to say special thanks to Péter Horváth, PhD (HA5CQA) and János Selmeczi, PhD (HA5FT) for their continous help and support.

Usage by example
----------------

### Demodulate WFM

    rtl_sdr -s 240000 -f 89500000 -g 20 - | csdr convert_u8_f | csdr fmdemod_quadri_cf | csdr fractional_decimator_ff 5 | csdr deemphasis_wfm_ff 48000 50e-6 | csdr convert_f_s16 | mplayer -cache 1024 -quiet -rawaudio samplesize=2:channels=1:rate=48000 -demuxer rawaudio -

- Baseband I/Q signal is coming from an RTL-SDR USB dongle, with a center frequency of `-f 104300000` Hz, a sampling rate of `-s 240000` samples per second.
- The `rtl_sdr` tool outputs an unsigned 8-bit I/Q signal (one byte of I sample and one byte of Q coming after each other), but `libcsdr` DSP routines internally use floating point data type, so we convert the data stream of `unsigned char` to `float` by `csdr convert_u8_f`.
- We want to listen one radio station at the frequency `-f 89500000` Hz (89.5 MHz).
- No other radio station is within the sampled bandwidth, so we send the signal directly to the demodulator. (This is an easy, but not perfect solution as the anti-aliasing filter at RTL-SDR DDC is too short.)
- After FM demodulation we decimate the signal by a factor of 5 to match the rate of the audio card (240000 / 5 = 48000).
- A de-emphasis filter is used, because pre-emphasis is applied at the transmitter to compensate noise at higher frequencies. The time constant for de-emphasis for FM broadcasting in Europe is 50 microseconds (hence the `50e-6`).
- Also, `mplayer` cannot play floating point audio, so we convert our signal to a stream of 16-bit integers.  

### Demodulate WFM: advanced

    rtl_sdr -s 2400000 -f 89300000 -g 20 - | csdr convert_u8_f | csdr shift_addition_cc -0.085 | csdr fir_decimate_cc 10 0.05 HAMMING | csdr fmdemod_quadri_cf | csdr fractional_decimator_ff 5 | csdr deemphasis_wfm_ff 48000 50e-6 | csdr convert_f_s16 | mplayer -cache 1024 -quiet -rawaudio samplesize=2:channels=1:rate=48000 -demuxer rawaudio -

- We want to listen to one radio station, but input signal contains multiple stations, and its bandwidth is too large for sending it directly to the FM demodulator.
- We shift the signal to the center frequency of the station we want to receive: `-0.085*2400000 = -204000`, so basically we will listen to the radio station centered at 89504000 Hz.
- We decimate the signal by a factor of 10. The transition bandwidth of the FIR filter used for decimation will be 10% of total bandwidth (as of parameter 0.05 is 10% of 0.5). Hamming window will be used for windowed FIR filter design.

Sample rates look like this:


                 2.4 Msps                     240 ksps                                  48 ksps
    I/Q source ------------> FIR decimation ------------> FM demod -> frac. decimation ---------> deemphasis -> sound card


*Note:* there is an example shell script that does this for you (without the unnecessary shift operation). If you just want to listen to FM radio, type:

    csdr-fm 89.5 20

The first parameter is the frequency in MHz, and the second optional parameter is the RTL-SDR tuner gain in dB.

### Demodulate NFM

    rtl_sdr -s 2400000 -f 145000000 -g 20 - | csdr convert_u8_f | csdr shift_addition_cc `python -c "print float(145000000-145350000)/2400000"` | csdr fir_decimate_cc 50 0.005 HAMMING | csdr fmdemod_quadri_cf | csdr limit_ff | csdr deemphasis_nfm_ff 48000 | csdr fastagc_ff | csdr convert_f_s16 | mplayer -cache 1024 -quiet -rawaudio samplesize=2:channels=1:rate=48000 -demuxer rawaudio -

- Note that the decimation factor is higher (we want to select a ~25 kHz channel).
- Also there is a python hack to calculate the relative shift offset. The real receiver frequency is `145350000` Hz.
- The de-emphasis filter is a fixed FIR filter that has a passband of 400-4000 Hz, also with a roll-off of -20 dB/decade.

### Demodulate AM

    rtl_sdr -s 2400000 -f 145000000 -g 20 - | csdr convert_u8_f | csdr shift_addition_cc `python -c "print float(145000000-144400000)/2400000"` | csdr fir_decimate_cc 50 0.005 HAMMING | csdr amdemod_cf | csdr fastdcblock_ff | csdr agc_ff | csdr limit_ff | csdr convert_f_s16 | mplayer -cache 1024 -quiet -rawaudio samplesize=2:channels=1:rate=48000 -demuxer rawaudio -

- `amdemod_cf` is used as demodulator.
- `agc_ff` should be used for AM and SSB.

### Design FIR band-pass filter (with complex taps)

    csdr firdes_bandpass_c 0 0.5 59 HAMMING --octave | octave -i

- ...and then plot its frequency response with octave. (You can close octave window by issuing Ctrl-C in the terminal window.)
- It will design a filter that lets only the positive frequencies pass (low cut is 0, high cut is 0.5 - these are relative to the sampling rate).
- If `--octave` and everything that follows is removed from the command, you get only the taps. E. g. the raw output of `firdes_lowpass_f` can be easily copied to C code.

### Demodulate SSB

    rtl_sdr -s 2400000 -f 145000000 -g 20 - | csdr convert_u8_f | csdr shift_addition_cc `python -c "print float(145000000-144400000)/2400000"` | csdr fir_decimate_cc 50 0.005 HAMMING | csdr bandpass_fir_fft_cc 0 0.1 0.05 | csdr realpart_cf | csdr agc_ff | csdr limit_ff | csdr convert_f_s16 | mplayer -cache 1024 -quiet -rawaudio samplesize=2:channels=1:rate=48000 -demuxer rawaudio -

- It is a modified Weaver-demodulator. The complex FIR filter removes the lower sideband and lets only the upper pass (USB). If you want to demodulate LSB, change `bandpass_fir_fft_cc 0 0.1` to `bandpass_fir_fft_cc -0.1 0`.

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
- `s16` is `signed short` of 2 bytes/16 bits (e. g. sound card input is usually `s16`)

Functions usually end as:

- `_ff` float input, float output
- `_cf` complex input, float output
- `_cc` complex input, complex output

Regarding *csdr*, it can convert a real/complex stream from one data format to another, to interface it with other SDR tools and the sound card.
The following commands are available:

- `csdr convert_u8_f`
- `csdr convert_f_u8`
- `csdr convert_s8_f`
- `csdr convert_f_s8`
- `csdr convert_s16_f`
- `csdr convert_f_s16`
- `csdr convert_s24_f [--bigendian]`
- `csdr convert_f_s24 [--bigendian]`

How to interpret: `csdr convert_<src>_<dst>`
You can use these commands on complex streams, too, as they are only interleaved values (I,Q,I,Q,I,Q... coming after each other).

> Note: The the functions with `i16` in their names have been renamed, but still work (e.g. `csdr convert_f_i16`).


### csdr commands

`csdr` should be considered as a reference implementation on using `libcsdr`. For additional details on how to use the library, check `csdr.c` and `libcsdr.c`.

Regarding `csdr`, the first command-line parameter is the name of a function, others are the parameters for the given function. Compulsory parameters are noted as `<parameter>`, optional parameters are noted as `[parameter]`.
Optional parameters have safe defaults, for more info look at the code.

----

### [realpart_cf](#realpart_cf)

Syntax:

    csdr realpart_cf

It takes the real part of the complex signal, and throws away the imaginary part.

----

### [clipdetect_ff](#clipdetect_ff)

Syntax:

    csdr clipdetect_ff

It clones the signal (the input and the output is the same), but it prints a warning on `stderr` if any sample value is out of the -1.0 ... 1.0 range.

----

### [limit_ff](#limit_ff)

Syntax:

    csdr limit_ff [max_amplitude]

The input signal amplitude will not be let out of the `-max_amplitude ... max_amplitude` range.

----

### [gain_ff](#gain_ff)

Syntax:

    csdr gain_ff <gain>

It multiplies all samples by `gain`.

----

### [clone](#clone)

Syntax:

    csdr clone

It copies the input to the output.

----

### [through](#through)

Syntax:

    csdr through

It copies the input to the output, while also displaying the data rate going through it.

----

### [none](#none)

Syntax:

    csdr none

The `csdr` process just exits with 0.

----

### [yes_f](#yes_f)

Syntax:

    csdr yes_f <to_repeat> [buf_times]

It outputs continously the `to_repeat` float number.

If `buf_times` is not given, it never stops.

Else, after outputing `buf_times` number of buffers (the size of which is stated in the `BUFSIZE` macro), it exits.

----

### [detect_nan_ff](#detect_nan_ff)

Syntax:

    csdr detect_nan_ff

Along with copying its input samples to the output, it prints a warning message to *stderr* if it finds any IEEE floating point NaN values among the samples.

----

### [dump_f](#dump_f)

Syntax:

    csdr dump_f

It prints all floating point input samples as text.

The C format string used is `"%g "`.

You can also use it to print complex float values, then you will see the I and Q samples interleaved, like: `I Q I Q I Q ...`

Alternatively, you can use the `od` command (built into most Linux distributions). To display a list of floating point values with their addresses as well, you can use: `od -vf`

Aliases for this function: `floatdump_f`

----

### [dump_u8](#dump_u8)

Syntax:

    csdr dump_u8

It prints all input bytes as text, in hexadecimal form. 

Alternatively, you can use the `xxd` command (built into most Linux distributions). To display a hexadecimal dump of the standard input (with addresses at the beginning of rows), you can use: `xxd -g1`

----

### [flowcontrol](#flowcontrol)

Syntax:

    csdr flowcontrol <data_rate> <reads_per_second>

It limits the data rate of a stream to a given `data_rate` number of bytes per second.

It copies `data_rate / reads_per_second` bytes from the input to the output, doing it `reads_per_second` times every second.

----

### [shift_math_cc](#shift_math_cc)

Syntax:

    csdr shift_math_cc <rate>

It shifts the signal in the frequency domain by `rate`.

`rate` is a floating point number between -0.5 and 0.5.

`rate` is relative to the sampling rate.

Internally, a sine and cosine wave is generated, and this function uses `math.h` for this purpose, which is quite accurate, but not always very fast.

----

### [shift_addition_cc](#shift_addition_cc)

Syntax:

    csdr shift_addition_cc <rate>

Operation is the same as for `shift_math_cc`.

Internally, this function uses trigonometric addition formulas to generate sine and cosine, which is a bit faster. (About 4 times on the machine I have tested it on.)

----

### [shift_addition_cc_test](#shift_addition_cc_test)

Syntax: 

    csdr shift_addition_cc_test

This function was used to test the accuracy of the method above.

----

### [shift_table_cc](#shift_table_cc)

Syntax: 

    csdr shift_table_cc <rate> [table_size]

Operation is the same as with `shift_math_cc`.

Internally, this function uses a look-up table (LUT) to recall the values of the sine function (for the first quadrant).

The higher the table size is, the smaller the phase error is.

----

### [shift_addfast_cc](#shift_addfast_cc)

Syntax: 

    csdr shift_addfast_cc <rate>

Operation is the same as for `shift_math_cc`.

Internally, this function uses a NEON-accelerated algorithm on capable systems, so it is advised to use this one on ARM boards.

----

### [shift_unroll_cc](#shift_unroll_cc)

Syntax: 

    csdr shift_unroll_cc <rate>

Operation is the same as for `shift_math_cc`.

This uses a modified algoritm that first stores a vector of sine and cosine values for given phase differences.

The loop in this function unrolls quite well if compiled on a PC. It was the fastest one on an i7 CPU during the tests.

----

### [decimating_shift_addition_cc](#decimating_shift_addition_cc)

Syntax: 

    csdr decimating_shift_addition_cc <rate> [decimation]

It shifts the input signal in the frequency domain, and also decimates it, without filtering. It will be useful as a part of the FFT channelizer implementation (to be done).

It cannot be used as a channelizer by itself, use `fir_decimate_cc` instead.

----

### [shift_addition_fc](#shift_addition_fc)

Syntax:

    csdr shift_addition_fc <rate>

It converts the real input signal to complex, and then shifts it in the frequency domain by `rate`.

----

### [dcblock_ff](#dcblock_ff)

Syntax: 

    csdr dcblock_ff

This is a DC blocking IIR filter.

----

### [fastdcblock_ff](#fastdcblock_ff)

Syntax: 

    csdr fastdcblock_ff

This is a DC blocker that works based on the average of the buffer.

----

### [fmdemod_atan_cf](#fmdemod_atan_cf)

Syntax: 

    csdr fmdemod_atan_cf

It is an FM demodulator that internally uses the `atan` function in `math.h`, so it is not so fast.

----

### [fmdemod_quadri_cf](#fmdemod_quadri_cf)

Syntax: 

    csdr fmdemod_quadri_cf

It is an FM demodulator that is based on the quadri-correlator method, and it can be effectively auto-vectorized, so it should be faster.

----

### [fmdemod_quadri_novect_cf](#fmdemod_quadri_novect_cf)

Syntax: 

    csdr fmdemod_quadri_novect_cf

It has more easily understandable code than the previous one, but can't be auto-vectorized.

----

### [deemphasis_wfm_ff](#deemphasis_wfm_ff)

Syntax: 

    csdr deemphasis_wfm_ff <sample_rate> <tau>

It does de-emphasis with the given RC time constant `tau`.

Different parts of the world use different pre-emphasis filters for FM broadcasting.

In Europe, `tau` should be chosen as `50e-6`, and in the USA, `tau` should be `75e-6`.

----

### [deemphasis_nfm_ff](#deemphasis_nfm_ff)

Syntax: 

    csdr deemphasis_nfm_ff <one_of_the_predefined_sample_rates>

It does de-emphasis on narrow-band FM for communication equipment (e.g. two-way radios).

It uses fixed filters so it works only on predefined sample rates, for the actual list of them run: 

    cat libcsdr.c | grep DNFMFF_ADD_ARRAY

----

### [amdemod_cf](#amdemod_cf)

Syntax: 

    csdr amdemod_cf

It is an AM demodulator that uses `sqrt`. On some architectures `sqrt` can be directly calculated by dedicated CPU instructions, but on others it may be slower.

----

### [amdemod_estimator_cf](#amdemod_estimator_cf)

Syntax: 

    csdr amdemod_estimator_cf

It is an AM demodulator that uses an estimation method that is faster but less accurate than `amdemod_cf`.

----

### [firdes_lowpass_f](#firdes_lowpass_f)

Syntax: 

    csdr firdes_lowpass_f <cutoff_rate> <length> [window [--octave]]

Low-pass FIR filter design function to output real taps, with a `cutoff_rate` proportional to the sampling frequency, using the windowed sinc filter design method.

`cutoff_rate` can be between 0 and 0.5.

`length` is the number of filter taps to output, and should be odd.

The longer the filter kernel is, the shorter the transition bandwidth is, but the more CPU time it takes to process the filter.

The transition bandwidth (proportional to the sampling rate) can be calculated as: `transition_bw = 4 / length`.

Some functions (below) require the `transition_bw` to be given instead of filter `length`. Try to find the best compromise between speed and accuracy by changing this parameter.

`window` is the window function used to compensate finite filter length. Its typical values are: `HAMMING`, `BLACKMAN`, `BOXCAR`. For the actual list of values, run: `cpp libcsdr.c | grep window\ ==`

The `--octave` parameter lets you directly view the filter response in `octave`. For more information, look at the [Usage by example] section.

----

### [firdes_bandpass_c](#firdes_bandpass_c)

Syntax: 

    csdr firdes_bandpass_c <low_cut> <high_cut> <length> [window [--octave]]

Band-pass FIR filter design function to output complex taps.

`low_cut` and ` high_cut` both may be between -0.5 and 0.5, and are also proportional to the sampling frequency.

Other parameters were explained above at `firdes_lowpass_f`.

----

### [fir_decimate_cc](#fir_decimate_cc)

Syntax: 

    csdr fir_decimate_cc <decimation_factor> [transition_bw [window]]

It is a decimator that keeps one sample out of `decimation_factor` samples.

To avoid aliasing, it runs a filter on the signal and removes spectral components above `0.5 × nyquist_frequency × decimation_factor` from the input signal.

----

### [fir_interpolate_cc](#fir_interpolate_cc)

Syntax: 

    csdr fir_interpolate_cc <interpolation_factor> [transition_bw [window]]

It is an interpolator that generates `interpolation_factor` number of output samples from one input sample.

To avoid aliasing, it runs a filter on the signal and removes spectral components above `0.5 × nyquist_frequency / interpolation_factor` from the output signal.

`transition_bw` and `window` are the parameters of the filter.

----

### [rational_resampler_ff](#rational_resampler_ff)

Syntax: 

    csdr rational_resampler_ff <interpolation> <decimation> [transition_bw [window]]

It is a resampler that takes integer values of `interpolation` and `decimation`.
The output sample rate will be `interpolation / decimation × input_sample_rate`.

`transition_bw` and `window` are the parameters of the filter.

----

### [fractional_decimator_ff](#fractional_decimator_ff)

Syntax: 

    csdr fractional_decimator_ff <decimation_rate> [num_poly_points ( [transition_bw [window]] | --prefilter )]

It can decimate by a floating point ratio.

It uses Lagrance interpolation, where `num_poly_points` (12 by default) input samples are taken into consideration while calculating one output sample. 

It can filter the signal with an anti-aliasing FIR filter before applying the Lagrange interpolation. This filter is inactive by default, but can be activated by:

* passing only the `transition_bw`, or both the `transition_bw` and the `window` parameters of the filter,
* using the `--prefilter` switch after `num_poly_points` to switch this filter on with the default parameters.

----

### [old_fractional_decimator_ff](#old_fractional_decimator_ff)

Syntax: 

    csdr old_fractional_decimator_ff <decimation_rate> [transition_bw [window]]

This is the deprecated, old algorithm to decimate by a floating point ratio, superseded by `fractional_decimator_ff`. 

(It uses linear interpolation, and its filter cuts at 59% of the passband.)

----

### [bandpass_fir_fft_cc](#bandpass_fir_fft_cc)

Syntax: 

    csdr bandpass_fir_fft_cc <low_cut> <high_cut> <transition_bw> [window]

It performs a bandpass FIR filter on complex samples, using FFT and the overlap-add method.

Parameters are described under `firdes_bandpass_c` and `firdes_lowpass_f`.

----

### [agc_ff](#agc_ff)

Syntax: 

    csdr agc_ff [hang_time [reference [attack_rate [decay_rate [max_gain [attack_wait [filter_alpha]]]]]]]

It is an automatic gain control function.

- `hang_time` is the number of samples to wait before starting to increase the gain after a peak.
- `reference` is the reference level for the AGC. It tries to keep the amplitude of the output signal close to that.
- `attack_rate` is the rate of decreasing the signal level if it gets higher than it used to be before.
- `decay_rate` is the rate of increasing the signal level if it gets lower than it used to be before.
- AGC won't increase the gain over `max_gain`.
- `attack_wait` is the number of sampels to wait before starting to decrease the gain, because sometimes very short peaks happen, and we don't want them to spoil the reception by substantially decreasing the gain of the AGC.
- `filter_alpha` is the parameter of the loop filter.

Its default parameters work best for an audio signal sampled at 48000 Hz.

----

### [fastagc_ff](#fastagc_ff)

Syntax: 

    csdr fastagc_ff [block_size [reference]]

It is a faster AGC that linearly changes the gain, taking the highest amplitude peak in the buffer into consideration. Its output will never exceed `-reference ... reference`.

----

### [fft_cc](#fft_cc)

Syntax: 

    csdr fft_cc <fft_size> <out_of_every_n_samples> [window [--octave] [--benchmark]]

It performs an FFT on the first `fft_size` samples out of `out_of_every_n_samples`, thus skipping `out_of_every_n_samples - fft_size` samples in the input.

It can draw the spectrum by using `--octave`, for more information, look at the [Usage by example] section.

FFTW can be faster if we let it optimalize a while before starting the first transform, hence the `--benchmark` switch.

----

### [fft_fc](#fft_fc)

Syntax:

    csdr fft_fc <fft_out_size> <out_of_every_n_samples> [--benchmark]

It works similarly to <a href="#fft_cc">fft_cc</a>, but on real input samples. 

For real FFT, the `fft_out_size` parameter is the number of output complex bins instead of the actual FFT size.

Number of input samples used for each FFT is `2 × fft_out_size`. This makes it easier to replace `fft_cc` by `fft_fc` in some applications.

----

### [fft_benchmark](#fft_benchmark)

Syntax: 

    csdr fft_benchmark <fft_size> <fft_cycles> [--benchmark]

It measures the time taken to process `fft_cycles` transforms of `fft_size`.
It lets FFTW optimalize if used with the `--benchmark` switch.

----

### [logpower_cf](#logpower_cf)

Syntax: 

    csdr logpower_cf [add_db]

Calculates `10*log10(i^2+q^2)+add_db` for the input complex samples. It is useful for drawing power spectrum graphs.

----

### [encode_ima_adpcm_i16_u8](#encode_ima_adpcm_i16_u8)

Syntax: 

    csdr csdr encode_ima_adpcm_i16_u8

Encodes the audio stream to IMA ADPCM, which decreases the size to 25% of the original.

----

### [decode_ima_adpcm_u8_i16](#decode_ima_adpcm_u8_i16)

Syntax: 

    csdr decode_ima_adpcm_u8_i16

Decodes the audio stream from IMA ADPCM.

----

### [compress_fft_adpcm_f_u8](#compress_fft_adpcm_f_u8)

Syntax: 

    csdr compress_fft_adpcm_f_u8 <fft_size>

Encodes the FFT output vectors of `fft_size`. It should be used on the data output from `logpower_cf`.

It resets the ADPCM encoder at the beginning of every vector, and to compensate it, `COMPRESS_FFT_PAD_N` samples are added at beginning (these equal to the first relevant sample).

The actual number of padding samples can be determined by running:

    cat csdr.c | grep "define COMPRESS_FFT_PAD_N"

----

### [fft_exchange_sides_ff](#fft_exchange_sides_ff)

Syntax: 

    csdr fft_exchange_sides_ff <fft_size>

It exchanges the first and second part of the FFT vector, to prepare it for the waterfall/spectrum display. It should operate on the data output from `logpower_cf`.

----

### [dsb_fc](#dsb_fc)

Syntax: 

    csdr dsb_fc [q_value]

It converts a real signal to a double sideband complex signal centered around DC.

It does so by generating a complex signal:
* the real part of which is the input real signal,
* the imaginary part of which is `q_value` (0 by default).

With `q_value = 0` it is an AM-DSB/SC modulator. If you want to get an AM-DSB signal, you will have to add a carrier to it.

----

### [add_dcoffset_cc](#add_dcoffset_cc)

Syntax: 

    csdr add_dcoffset_cc

It adds a DC offset to the complex signal: `i_output = 0.5 + i_input / 2, q_output = q_input / 2`

----

### [convert_f_samplerf](#convert_f_samplerf)

Syntax: 

    csdr convert_f_samplerf <wait_for_this_sample>

It converts a real signal to the `-mRF` input format of [https://github.com/F5OEO/rpitx](rpitx), so it allows you to generate frequency modulation. The input signal will be the modulating signal. The `<wait_for_this_sample>` parameter is the value for `rpitx` indicating the time to wait between samples. For a sampling rate of 48 ksps, this is 20833.

----

### [fmmod_fc](#fmmod_fc)

Syntax: 

    csdr fmmod_fc

It generates a complex FM modulated output from a real input signal.

----

### [fixed_amplitude_cc](#fixed_amplitude_cc)

Syntax: 

    csdr fixed_amplitude_cc <new_amplitude>

It changes the amplitude of every complex input sample to a fixed value. It does not change the phase information of the samples.

----

### [mono2stereo_s16](#mono2stereo_s16)

Syntax: 

    csdr mono2stereo_s16

It doubles every input sample. 

Example: if the input samples are 16 bit signed integers: 

    23 -452 3112

The output will be:

    23 23 -452 -452 3112 3112

----

### [setbuf](#setbuf)

Syntax: 

    csdr setbuf <buffer_size>

See the [buffer sizes](#buffer_sizes) section.

    squelch_and_smeter_cc --fifo <squelch_fifo> --outfifo <smeter_fifo> <use_every_nth> <report_every_nth>

This is a controllable squelch, which reads the squelch level input from `<squelch_fifo>` and writes the power level output to `<smeter_fifo>`. Both input and output are in the format of `%g\n`. While calculating the power level, it takes only every `<use_every_nth>` sample into consideration. It writes the S-meter value for every `<report_every_nth>` buffer to `<smeter_fifo>`. If the squelch level is set to 0, it it forces the squelch to be open. If the squelch is closed, it fills the output with zero.

----

### [fifo](#fifo)

Syntax: 

    csdr fifo <buffer_size> <number_of_buffers>

It is similar to `clone`, but internally it uses a circular buffer. It reads as much as possible from the input. It discards input samples if the input buffer is full.

----

### [psk31_varicode_encoder_u8_u8](#psk31_varicode_encoder_u8_u8)

Syntax: 

    csdr psk31_varicode_encoder_u8_u8

It encodes ASCII characters into varicode for PSK31 transmission. It puts a `00` sequence between the varicode characters (which acts as a separator). 

For the Varicode character set, see: http://www.arrl.org/psk31-spec

For example, `aaa` means the bit sequence `101100101100101100`. 

For this input, the output of `psk31_varicode_encoder_u8_u8` will be the following bytes (in hexadecimal):

```
01 00 01 01 00 00 01 00 
01 01 00 00 01 00 01 01
00 00
```

----

### [repeat_u8](#repeat_u8)

Syntax:

    csdr repeat_u8 <data_bytes × N>

It repeatedly outputs a set of data bytes (given with decimal numbers).

For example, `csdr repeat_u8 1 1 0 0` will output:

```
01 01 00 00 01 01 00 00 
01 01 00 00 01 01 00 00
```

----

### [uniform_noise_f](#uniform_noise_f)

Syntax:

    csdr uniform_noise_f

It outputs uniform white noise. All samples are within the range [-1.0, 1.0].

----

### [gaussian_noise_c](#gaussian_noise_c)

Syntax:

    csdr gaussian_noise_c

It outputs Gaussian white noise. All samples are within the unit circle.

----

### [pack_bits_8to1_u8_u8](#pack_bits_8to1_u8_u8)

Syntax:

    csdr pack_bits_8to1_u8_u8 

TODO 

----

### [pack_bits_1to8_u8_u8](#pack_bits_1to8_u8_u8)

Syntax:

    csdr pack_bits_1to8_u8_u8 

It serializes the bytes on the input: it outputs each bit of the input byte as a single byte valued 0x00 or 0x01, starting from the lowest bit and going to the highest bit.

The output is 8 times as large in size as the input. 

For example, the input byte 0x43 will result in eight bytes at the output:

```
01 01 00 00 00 00 01 00 
```

For consequtive 0x02, 0x03, 0xff bytes on the input, the output will be:

```
00 01 00 00 00 00 00 00 
01 01 00 00 00 00 00 00 
01 01 01 01 01 01 01 01 
```

----

### [awgn_cc](#awgn_cc)

Syntax:

    csdr awgn_cc <snr_db> [--snrshow]

It adds white noise with the given SNR to a signal assumed to be of 0 dB power.

If the `--snrshow` switch is given, it also shows the actual SNR based on the calculated power of signal and noise components.

----

### [add_n_zero_samples_at_beginning_f](#add_n_zero_samples_at_beginning_f)

Syntax:

    csdr add_n_zero_samples_at_beginning_f <n_zero_samples> 

When the function is executed, it furst writes `<n_zero_samples>` 32-bit floating point zeros at the output, after that it just clones the input at the output. 

----

### [fft_one_side_ff](#fft_one_side_ff)

Syntax:

    csdr fft_one_side_ff <fft_size>

If the frequency domain signal spans between frequencies -fs/2 to fs/2, this function removes the part from -fs/2 to DC. This can be useful if the FFT of a real signal has been taken (so that the spectrum is mirrored to DC).  

----

### [logaveragepower_cf](#logaveragepower_cf)

Syntax:

    csdr logaveragepower_cf <add_db> <fft_size> <avgnumber>

It works like <a href="#logpower_cf">logpower_cf </a>, but it calculates the average of every `avgnumber` FFTs. 

----

### [mono2stereo_s16](#mono2stereo_s16)

Syntax:

    csdr mono2stereo_s16

It duplicates each 16-bit integer input sample.

----

### [psk31_varicode_decoder_u8_u8](#psk31_varicode_decoder_u8_u8)

Syntax:

    csdr psk31_varicode_decoder_u8_u8

It expects symbols encoded as 0x00 and 0x01 bytes on the input, and extracts Varicode characters from them.  

----

### [_fft2octave](#_fft2octave)

Syntax:

    csdr _fft2octave <fft_size>

It is used for plotting FFT data with a GNU Octave session, piping its output to `octave -i`.

----

### [invert_u8_u8](#invert_u8_u8)

Syntax:

    csdr invert_u8_u8

It maps

* each 0x00 to 0x01,
* each 0x01 to 0x00.

----

### [rtty_baudot2ascii_u8_u8](#rtty_baudot2ascii_u8_u8)

Syntax:

    csdr rtty_baudot2ascii_u8_u8

This function awaits baudot code characters on its input (ranging from 0b00000000 to 0b00011111), and converts them into ASCII characters. It has an internal state to switch between letters and figures. 

----

### [binary_slicer_f_u8](#binary_slicer_f_u8)

Syntax:

    csdr binary_slicer_f_u8

* If the input sample is below or equals to 0.0, it outputs a 0x00.
* If the input sample is above 0.0, it outputs a 0x01.

----

### [serial_line_decoder_f_u8](#serial_line_decoder_f_u8)

Syntax:

    csdr serial_line_decoder_f_u8 <samples_per_bits> [databits [stopbits]]

It decodes bits from a sampled serial line. It does so by finding the appropriate start and stop bits, and extracts the data bits in between.

----

### [pll_cc](#pll_cc)

Syntax:

    csdr pll_cc (1 [alpha] |2 [bandwidth [damping_factor [ko [kd]]]])

It implements a PLL that can lock onto a sinusoidal input signal. 

The first parameter corresponds to the order of the PLL loop filter (first or second order), others are parameters of the loop filter.

----

### [timing_recovery_cc](#timing_recovery_cc)

Syntax:

    csdr timing_recovery_cc  <algorithm> <decimation> [mu [max_error [--add_q [--output_error | --output_indexes |  --octave <show_every_nth> |  --octave_save <show_every_nth> <directory> ]]]] 

It implements non-data aided timing recovery (Gardner and early-late gate algorithms). 

[More information](http://openwebrx.org/msc-thesis.pdf#page=34) (section 4.4 from page 34)

----

### [octave_complex_c](#octave_complex_c)

Syntax:

    csdr octave_complex_c <samples_to_plot> <out_of_n_samples> [--2d]

It generates octave commands to plot a complex time domain signal. Its output can be piped into `octave -i`. It plots every `samples_to_plot` samples `out_of_n_samples`.

----

### [psk_modulator_u8_c](#psk_modulator_u8_c)

Syntax:

    csdr psk_modulator_u8_c <n_psk>

It generates an N-PSK modulated signal from the input symbols. 

As an example, for `n_psk`=4, it will translate:

* any 0x00 byte on the input into 1+0j on the output,
* any 0x01 byte on the input into 0+1j on the output,
* any 0x02 byte on the input into -1+0j on the output,
* any 0x03 byte on the input into 0-1j on the output.

----

### [duplicate_samples_ntimes_u8_u8](#duplicate_samples_ntimes_u8_u8)

Syntax:

    csdr duplicate_samples_ntimes_u8_u8 <sample_size_bytes> <ntimes>

It duplicates each sample of `sample_size_bytes` the given `ntimes` times.

----

### [psk31_interpolate_sine_cc](#psk31_interpolate_sine_cc)

Syntax:

    csdr psk31_interpolate_sine_cc <interpolation>

The input to this function is one complex sample per symbol, the output is `interpolation` samples per symbol, interpolated using a cosine envelope (which is used for PSK31). 

----

### [differential_encoder_u8_u8](#differential_encoder_u8_u8)

Syntax:

    csdr differential_encoder_u8_u8

It can be used while generating e.g. differential BPSK modulation. 

* If the input is 0x01, the output remains the same as the last output.
* If the input is 0x00, the output changes from 0x00 to 0x01, or 0x01 to 0x00.

----

### [differential_decoder_u8_u8](#differential_decoder_u8_u8)

Syntax:

    csdr differential_decoder_u8_u8

It can be used while demodulating e.g. differential BPSK modulation. The following table show the logic function it performs:

| Last input | Current input | Output | 
| ---------- | ------------- | ------ |
| 0x00       | 0x00          | 0x01   |
| 0x00       | 0x01          | 0x00   |
| 0x01       | 0x00          | 0x00   |
| 0x01       | 0x01          | 0x01   |

----

### [bpsk_costas_loop_cc](#bpsk_costas_loop_cc)

Syntax:

    csdr bpsk_costas_loop_cc <loop_bandwidth> <damping_factor> [--dd | --decision_directed] [--output_error | --output_dphase | --output_nco | --output_combined <error_file> <dphase_file> <nco_file>]

It implements a Costas loop for BPSK signals. 

[More information](http://openwebrx.org/msc-thesis.pdf#page=55) (section 5.4 from page 55)

----

### [simple_agc_cc](#simple_agc_cc)

Syntax:

    csdr simple_agc_cc <rate> [reference [max_gain]] 

It is an automatic gain control function with a single pole IIR loop filter.

- `reference` is the reference level for the AGC. It tries to keep the amplitude of the output signal close to that.
- AGC won't increase the gain over `max_gain`.
- `rate` is the parameter of the loop filter.

The block diagram of this function follows:

![simple_agc_cc block diagram](https://raw.githubusercontent.com/wiki/simonyiszk/csdr/simple-agc-dataflow.png)

----

### [peaks_fir_cc](#peaks_fir_cc)

Syntax:

    csdr peaks_fir_cc <taps_length> <peak_rate × N>

It applies a peak filter to the input signal. The peak filter is a very narrow bandpass filter, the opposite of a notch filter. The higher the `taps_length` is, the sharper the filter frequency transfer function is. 

`peak_rate` is the center of the passband, in proportion to the sampling rate. 

----

### [firdes_peak_c](#firdes_peak_c)

Syntax:

    csdr firdes_peak_c <rate> <length> [window [--octave]]

It designs a FIR peak filter, and writes the taps to the output. More about this filter at <a href="#peaks_fir_cc">peaks_fir_cc</a>. 

This command also supports GNU Octave-friendly output that can be piped into the Octave interpreter `octave -i`.

----

### [normalized_timing_variance_u32_f](#normalized_timing_variance_u32_f)

Syntax:

    csdr normalized_timing_variance_u32_f <samples_per_symbol> <initial_sample_offset> [--debug]

It calculates the normalized timing variance. It works on the sample indexes output from the `timing_recovery_cc` function. 

----

### [pulse_shaping_filter_cc](#pulse_shaping_filter_cc)

Syntax:

    csdr pulse_shaping_filter_cc (RRC <samples_per_symbol> <num_taps> <beta> | COSINE <samples_per_symbol>)

It runs a pulse shaping FIR filter on the signal.

* `RRC` stands for Root-Raised-Cosine filter, a design parameter of which is `beta`.
* The `COSINE` filter is the one used for BPSK31. 
* `samples_per_symbol` is the number of input samples per symbol.
* `num_taps` is the filter length.

----

### [firdes_pulse_shaping_filter_f](#firdes_pulse_shaping_filter_f)

Syntax:

    csdr firdes_pulse_shaping_filter_f (RRC <samples_per_symbol> <num_taps> <beta> | COSINE <samples_per_symbol>)

It designs a pulse shaping filter, and outputs the taps. It has the same parameters as `pulse_shaping_filter_cc`. 

----

### [generic_slicer_f_u8](#generic_slicer_f_u8)

Syntax:

    csdr generic_slicer_f_u8 <n_symbols>

It decides which symbol the sample corresponds to, where the highest symbol corresponds to 1.0, and the lowest symbol corresponds to -1.0.

As an example, if N=3, the 3 symbols to choose from are: -1, 0, 1. The algorithm will output:

* 0x00 for any input sample between -infinity and -0.5.
* 0x01 for any input sample between -0.5 and 0.5.
* 0x02 for any input sample between 0.5 and infinity. 

----

### [plain_interpolate_cc](#plain_interpolate_cc)

Syntax:

    csdr plain_interpolate_cc <interpolation>

It interpolates the signal by writing `interpolation - 1` zero samples between each input sample. You need to run an anti-aliasing filter on its output.

----

### [dbpsk_decoder_c_u8](#dbpsk_decoder_c_u8)

Syntax:

    csdr dbpsk_decoder_c_u8

It implements a differential BPSK demodulator, with the following data flow:

![DBPSK dataflow](https://raw.githubusercontent.com/wiki/simonyiszk/csdr/dbpsk-dataflow.png)

The output is 0x00 or 0x01.

----

### [bfsk_demod_cf](#bfsk_demod_cf)

Syntax:

    csdr bfsk_demod_cf <spacing> <filter_length>

It implements a 2-FSK demodulator, with the following data flow:

![BFSK dataflow](https://raw.githubusercontent.com/wiki/simonyiszk/csdr/bfsk-dataflow.png)

You can calculate the expected frequencies of the two tones on the input by the following formulas: `+(spacing/sampling_rate)` and `-(spacing/sampling_rate)`.

Filter length is the length of the peak filters (FIR) applied to the input för each tone.

----

### [add_const_cc](#add_const_cc)

Syntax:

    csdr add_const_cc <i> <q>

It adds a constant value of `i+q*j` to each input sample.

----

### [pattern_search_u8_u8](#pattern_search_u8_u8)

Syntax: 

    csdr pattern_search_u8_u8 <values_after> <pattern_values × N>

It can be used for preamble search. It looks for a given sequence of N bytes (`<pattern_values × N>`) in the input data, and if the sequence is found, it reads the following `<values_after>` bytes and outputs them. The `<pattern_values × N>` parameter is read as unsigned integers.

----

### [tee](#tee)

Syntax:

    csdr tee <path> [buffers]

Similarly to the `tee` command, it reads data from the standard input, and writes it to both a file and the standard output. 

Unlike `tee`, if it fails to flush the data to the file, it still flushes it to the standard output. This allows us to have less glitches / better response time if we use this as a way to put branches in the data flow. Example:

    mkfifo /tmp/csdr_fifo
    rtl_sdr - | csdr tee /tmp/csdr_fifo | csdr dump_u8
    cat /tmp/csdr_fifo | csdr convert_u8_f | csdr dump_f

How the data flow looks like:

    rtl_sdr --> tee --> dump_u8
                 |
                \/
             convert_u8_f --> dump_f


### [?](#search_the_function_list)

Syntax: 

    csdr ?<search_the_function_list>

You can search the functions available in `csdr` just as if you typed: `csdr 2>&1 | grep <search_what>`

### [=](#evaluate-python-expression)

Syntax:

    csdr =<evaluate_python_expression>

When running complicated `csdr` commands, we usually run into using `python` to calculate certain parameters. 

This function can eliminate some typing and make our command clearer.

Instead of having to write: 

    csdr shift_addition_cc $(python -c "print 1200/2400000.") 

...we can type: 

    csdr shift_addition_cc $(csdr =1200/2400000.)

If using parenthesis inside the expression, it needs to be escaped (as `bash` would want to parse it): 

    csdr shift_addition_cc $(csdr =\(1200+300\)/2400000)

Another solution is using single quotes to wrap the expression: 

    csdr shift_addition_cc $(csdr '=(1200+300)/2400000.')

Current version of `csdr` executes the following python script for this function:

```python
import os, sys
from math import *
print <evaluate_python_expression>
```

This means that one can also call math functions like `sqrt()`.

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

#### Buffer sizes

*csdr* has three modes of determining the buffer sizes, which can be chosen by the appropriate environment variables:
* *default:* 16k or 1k buffer is chosen based on function,
* *dynamic buffer size determination:* input buffer size is recommended by the previous process, output buffer size is determined by the process,
* *fixed buffer sizes*.

*csdr* can choose from two different buffer sizes by **default**.
* For operations handling the full-bandwidth I/Q data from the receiver, a buffer size of 16384 samples is used (see `env_csdr_fixed_big_bufsize` in the code).
* For operations handling only a selected channel, a buffer size of 1024 samples is used (see `env_csdr_fixed_bufsize` in the code).

*csdr* now has an experimental feature called **dynamic buffer size determination**, which is switched on by issuing `export CSDR_DYNAMIC_BUFSIZE_ON=1` in the shell before running `csdr`. If it is enabled:
* All `csdr` processes in a DSP chain acquire their recommended input buffer size from the previous `csdr` process. This information is in the first 8 bytes of the input stream.
* Each process can decide whether to use this or choose another input buffer size (if that's more practical).
* Every process sends out its output buffer size to the next process. Then it startss processing data.
* The DSP chain should start with a `csdr setbuf <buffer_size>` process, which only copies data from the input to the output, but also sends out the given buffer size information to the next process.
* The 8 bytes of information included in the beginning of the stream is:
  * a preamble of the bytes 'c','s','d','r' (4 bytes),
  * the buffer size stored as `int` (4 bytes).
* This size always counts as samples, as we expect that the user takes care of connecting the functions with right data types to each other.

> I added this feature while researching how to decrease the latency of a DSP chain consisting of several multirate algorithms.<br />
> For example, a `csdr fir_decimate_cc 10` would use an input buffer of 10240, and an output buffer of 1024. The next process in the chain, `csdr bandpass_fir_fft_cc` would automatically adjust to it, using a buffer of 1024 for both input and output.<br />
> In contrast to original expectations, using dynamic buffer sizes didn't decrease the latency much.

If dynamic buffer size determination is disabled, you can still set a **fixed buffer size** with `export CSDR_FIXED_BUFSIZE=<buffer_size>`.

For debug purposes, buffer sizes of all processes can be printed using `export CSDR_PRINT_BUFSIZES=1`.

If you add your own functions to `csdr`, you have to initialize the buffers before doing the processing. Buffer size will be stored in the global variable `the_bufsize`.

Example of initialization if the process generates N output samples for N input samples:

    if(!sendbufsize(initialize_buffers())) return -2;

Example of initalization if the process generates N/D output samples for N input samples:

    if(!initialize_buffers()) return -2;
    sendbufsize(the_bufsize/D);

Example of initialization if the process allocates memory for itself, and it doesn't want to use the global buffers:

    getbufsize(); //dummy
    sendbufsize(my_own_bufsize);

Example of initialization if the process always works with a fixed output size, regardless of the input:

    if(!initialize_buffers()) return -2;
    sendbufsize(fft_size);

#### Testbench

`csdr` was tested with GNU Radio Companion flowgraphs. These flowgraphs are available under the directory `grc_tests`, and they require the <a href="https://github.com/simonyiszk/gr-ha5kfu">gr-ha5kfu</a> set of blocks for GNU Radio.  

## [sdr.js](#sdrjs)

*sdr.js* is *libcsdr* compiled to JavaScript code with *Emscripten*. Nowadays JavaScript runs quite fast in browsers, as all major browser vendors included JavaScript JIT machines into their product. You can find a <a href="https://kripken.github.io/mloc_emscripten_talk/cppcon.html">great introductory slideshow here</a> on the concept behind *Emscripten* and *asm.js*.

The purpose of *sdr.js* is to make SDR DSP processing available in the web browser. However, it is not easy to use in production yet. By now, only those functions have wrappers that the front-end of OpenWebRX uses.

To compile *sdr.js*, first get <a href="http://emscripten.org/">Emscripten</a>. (It turns out that there is an *emscripten* package in Ubuntu repositories.)

To install and build dependencies (for now, only FFTW3):

    make emcc-get-deps

To compile *sdr.js* (which will be created in the `sdr.js` subdirectory):

    make emcc

You can test *sdr.js* by opening *sdr.html*. It contains a test for *firdes_lowpass_f* for this time.

To remove *sdr.js* and the compiled dependencies:

    make emcc-clean

## [nmux](#nmux)

The repo also contains a command line tool called `nmux`, which is a TCP stream multiplexer. It reads data from the standard input, and sends it to each client connected through TCP sockets. Available command line options are:
* `--port (-p), --address (-a):` TCP port and address to listen.
* `--bufsize (-b), --bufcnt (-n)`: Internal buffer size and count.
* `--help (-h)`: Show help message.

`nmux` was originally written for use in OpenWebRX.

## [Licensing](#licensing)

Most of the code of `libcsdr` is under BSD license.  
However, before the implementation of some algoritms, GPL-licensed code from other applications have been reviewed.
In order to eliminate any licesing issues, these parts are placed under a different file.
However, the library is still fully functional with BSD-only code, altough having only less-optimized versions of some algorithms.  
It should also be noted that if you compile with `-DUSE_FFTW` and `-DLIBCSDR_GPL` (as default), the GPL license would apply on the whole result.
