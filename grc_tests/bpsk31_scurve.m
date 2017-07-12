#!/usr/bin/octave

%{
function [output]=fgc(path, type)
	if(type(1)=='f')
	elseif(type(1)=='c')
	end
end
%}

function output=shrunf(cmd)
	SIGTERM=15;
	output=[];
	[pin, pout, pid]=popen2('bash',{'-c', cmd});
	%fclose(pin);
	sleep(0.1)
	do
		current_output=fread(pout, Inf, 'float32');
		output=[output; current_output];
	until(feof(pout))
	waitpid(pid);
	%kill(pid, SIGTERM);
    fclose(pin);
    fclose(pout);
end

function error_value=run_tr(skip, which_ted) 
	out_vect=shrunf(sprintf('dd bs=8 skip=%d if=bpsk31_baseband_sample_complex_8000_sps_010101.raw | csdr timing_recovery_cc %s 256 --add_q --output_error', skip, which_ted));
	error_value=out_vect(2);
end

function error_values=mkscurve(which_ted, skips)
    error_values=[]
    for skip=skips
        error_values=[error_values run_tr(skip, which_ted)];
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

skips_gardner=0:16:256
error_values_gardner=mkscurve('GARDNER',skips_gardner);
skips_earlylate=0:16:256
error_values_earlylate=mkscurve('EARLYLATE',skips_earlylate);

%graphics_toolkit("gnuplot")
h=figure(1);

plot((skips_gardner-128)/256, -error_values_gardner, 'linewidth', 2);
title('S-curve for Gardner TED');
fmtplot(h)
grid on
pause

plot((skips_earlylate-128)/256, error_values_earlylate, 'linewidth', 2);
title('S-curve for early-late TED');
fmtplot(h)
grid on
pause
