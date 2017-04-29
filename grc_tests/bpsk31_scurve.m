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
end

function error_value=run_tr(skip)
	out_vect=shrunf(sprintf('dd bs=8 skip=%d if=bpsk31_baseband_sample_complex_8000_sps_010101.raw | csdr timing_recovery_cc EARLYLATE 256 --add_q --output_error',skip));
	error_value=out_vect(2);
end

skips=0:256;
error_values=[]
for skip=skips
	error_values=[error_values run_tr(skip)];
end

error_values

%graphics_toolkit("gnuplot")
h=figure(1);
plot(skips, error_values, 'linewidth', 2);
xlabel('Phase offset in number of samples');
ylabel('Error value (TED output)');

FN = findall(h,'-property','FontName');
set(FN,'FontName','/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerifCondensed.ttf');
set(FN,'FontName','times');
FS = findall(h,'-property','FontSize');
set(FS,'FontSize',18);

pause
