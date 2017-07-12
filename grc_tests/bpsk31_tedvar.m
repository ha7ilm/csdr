#!/usr/bin/octave

%you need to first install the parallel and struct packages:
%pkg install -forge struct
%pkg install -forge parallel
pkg load parallel

function y=inarg(x)
    for i=1:length(argv())
        if strcmp(argv(){i},x)
            y=1;
            return
        end
    end
    y=0;
end

bpfcmd="csdr bandpass_fir_fft_cc $(csdr \"=-31.25/8e3\") $(csdr \"=31.25/8e3\") $(csdr \"=31.25/8e3\") | ";

if !inarg('--nogen')
    fwrite(stdout, "===========================================\nGenerating baseband signal from random data\n===========================================\n");
    system(["cat /dev/urandom | csdr pack_bits_8to1_u8_u8 | csdr psk_modulator_u8_c 2 | csdr gain_ff 0.25 | csdr psk31_interpolate_sine_cc 256 | " bpfcmd "csdr add_n_zero_samples_at_beginning_f 170 | pv -ps 2g | dd iflag=fullblock bs=128M count=16 of=/tmp/psk31-raw-data"]);
    fwrite(stdout, "===========================================\nGenerating Gaussian white noise for agwn_cc\n===========================================\n");
    system(["csdr gaussian_noise_c | " bpfcmd "pv -ps 256m | dd of=/tmp/psk31-gaussian-noise iflag=fullblock bs=256M count=1"]);
end
if inarg('--onlygen') 
    exit(0)
end
fwrite(stdout, "===========================================\nCalculating variance graph data            \n===========================================\n");

function output=shrun(cmd, type, minsize)
    SIGTERM=15;
    output=[];
    cmd
    [pin, pout, pid]=popen2('bash',{'-c', cmd});
    %fclose(pin);
    do
        sleep(0.3)
        fwrite(stdout,'.');
        %size(output)
        %output
        current_output=fread(pout, Inf, type);
        frewind(pout);
        output=[output; current_output];
    until(size(output)(1)>=minsize)
    waitpid(pid);
    kill(pid, SIGTERM);
    fclose(pin);
    fclose(pout);
end

function variance=run_var(snr, which_ted)
    disp('ran a command')
    out_vect=shrun(sprintf('cat /tmp/psk31-raw-data | csdr awgn_cc %d --awgnfile /tmp/psk31-gaussian-noise | csdr simple_agc_cc 0.0001 0.5 | csdr timing_recovery_cc %s 256 0.5 2 --add_q --output_indexes | CSDR_FIXED_BUFSIZE=1048576 csdr normalized_timing_variance_u32_f 256 85', snr, which_ted), 'float32', 1);
    disp('run_var output:');
    out_vect'
    variance=out_vect(1);
end

function variances=mkvarplot(which_ted, snrs)
    fun = @(x) run_var(x, which_ted);
    variances=pararrayfun(nproc, fun, snrs);
    %{
    variances=[]
    for snr=snrs
        snr
        variances=[variances run_var(snr, which_ted)];
    end
    %}
end

function fmtplot(h)
    FN = findall(h,'-property','FontName');
    set(FN,'FontName','/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerifCondensed.ttf');
    set(FN,'FontName','times');
    FS = findall(h,'-property','FontSize');
    set(FS,'FontSize',18);
    xlabel('E_b/N_0 [dB]');
    ylabel('Phase error variance [rad^2]');
end

%snrs=-10:5:10
snrs=-10:5:25
%snrs=[10]
error_values=mkvarplot('EARLYLATE',snrs);

%graphics_toolkit("gnuplot")
h=figure(1);

ebn0=snrs+9.7

semilogy(ebn0, error_values, 'linewidth', 2);
title('Estimation variance');
fmtplot(h)
pause

if !inarg('--nogen')
    system('rm /tmp/psk31-raw-data /tmp/psk31-gaussian-noise');
end
