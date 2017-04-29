#!/usr/bin/octave

%you need to first install the parallel and struct packages:
%pkg install -forge struct
%pkg install -forge parallel
pkg load parallel

system('cat /dev/urandom | csdr pack_bits_8to1_u8_u8 | csdr psk_modulator_u8_c 2 | csdr gain_ff 0.25 | csdr psk31_interpolate_sine_cc 256 | csdr add_n_zero_samples_at_beginning_f 170 | dd bs=32M count=1 of=/tmp/psk31-raw-data');

function output=shrun(cmd, type, minsize)
    SIGTERM=15;
    output=[];
    cmd
    [pin, pout, pid]=popen2('bash',{'-c', cmd});
    %fclose(pin);
    do
        sleep(0.3)
        disp('size(output)');
        %size(output)
        %output
        current_output=fread(pout, Inf, type)
        frewind(pout)
        output=[output; current_output];
    until(size(output)(1)>=minsize)
    waitpid(pid);
    kill(pid, SIGTERM);
    fclose(pin);
    fclose(pout);
end

function variance=run_var(snr, which_ted)
    disp('ran a command')
    out_vect=shrun(sprintf('cat /tmp/psk31-raw-data | csdr awgn_cc %d | csdr timing_recovery_cc %s 256 --add_q --output_indexes | CSDR_FIXED_BUFSIZE=65536 csdr normalized_timing_variance_u32_f 256 85', snr, which_ted), 'float32', 1);
    disp('run_var output:');
    out_vect'
    variance=out_vect(1);
end

function variances=mkvarplot(which_ted, snrs)
    %{
    fun = @(x) run_var(x, which_ted);
    variances=pararrayfun(nproc, fun, snrs);
    %}
    variances=[]
    for snr=snrs
        snr
        variances=[variances run_var(snr, which_ted)];
    end
end

function fmtplot(h)
    FN = findall(h,'-property','FontName');
    set(FN,'FontName','/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerifCondensed.ttf');
    set(FN,'FontName','times');
    FS = findall(h,'-property','FontSize');
    set(FS,'FontSize',18);
    xlabel('Phase offset in number of samples');
    ylabel('Error value (TED output)');
end

snrs_gardner=-30:5:40
error_values_gardner=mkvarplot('GARDNER',snrs_gardner);
%{
snrs_earlylate=0:256
error_values_earlylate=mkvarplot('EARLYLATE',snrs_earlylate);
%}

%graphics_toolkit("gnuplot")
h=figure(1);

semilogy(snrs_gardner, error_values_gardner, 'linewidth', 2);
title('S-curve for Gardner TED');
fmtplot(h)
pause

%{
semilogy(snrs_earlylate, error_values_earlylate, 'linewidth', 2);
title('S-curve for early-late TED');
fmtplot(h)
pause
%}

system('rm /tmp/psk31-raw-data');
