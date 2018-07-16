function phs=fractpitch2(fftl)
%	Phase rotator for fractional pitch
%	This program produces 'phs' as the phase rotator.

%	by Hideki Kawahara
%	22/August/1996

amp=15;
t=((1:fftl)-fftl/2-1)/fftl*2;
phs=t+(1-exp(amp*t))./(1+exp(amp*t)) ...
   -(1+(1-exp(amp))/(1+exp(amp)))*t;
phs(1)=0;
phs=phs*pi;