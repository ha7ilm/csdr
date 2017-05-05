#!/usr/bin/octave

global Tb=20

function g=gbb(t) %impulse response of pulse shaping filter
	global Tb
	g=t;
	for i = 1:size(t)(2)
		if (t(i)>1*Tb || t(i)<=-1*Tb)
			g(i) = 0;
		else
			g(i) = 0.5+cos((t(i)/(Tb*1))*pi)/2; %this is not RRC, rather a sinusoidal pulse shape
		end
	end
end

function [toreturny, plotrange]=y(s)
	global Tb
	slen=size(s)(2)
	plotrange=(-3*Tb):(slen+2)*Tb-1;
	plotlen=size(plotrange)(2)
	toreturny=zeros(1,plotlen);
	for i=1:slen %sum of (symbol[i] * filter impulse response) for all symbols
		toreturny+=s(i)*gbb(plotrange.-(i-1)*Tb); 
	end
end

function fmtplot(h)
    FN = findall(h,'-property','FontName');
    set(FN,'FontName','/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerifCondensed.ttf');
    set(FN,'FontName','times');
    FS = findall(h,'-property','FontSize');
    set(FS,'FontSize',18);
    set(FS,'FontSize',18);
end

h=figure(1);
subplot(2, 1, 1);
plot(y([1]), 'linewidth', 2)
title(sprintf("Impulse response of pulse shaping filter\n(1 symbol = %d samples)", Tb))
xlabel('t [1/Ts]')
ylabel('h(t)')

subplot(2, 1, 2);
plot(y([1 1 -1 -1 1 1 1 -1 1 -1 1 1]), 'linewidth', 2)
title("Baseband signal for modulator input\nbit sequence: 110011101011") %assuming that differential encoding has already been performed
xlabel('t [1/Ts]')
ylabel('s(t)')
fmtplot(h);
pause
