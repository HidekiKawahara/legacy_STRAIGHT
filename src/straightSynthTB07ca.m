function [sy,synthSataus]=straightSynthTB07ca(n2sgram,f0raw,shiftm,fs, ...
  pcnv,fconv,sconv,gdbw,delfrac,delsp,cornf,delfracind,ap,imap,imgi,lowestF0)
%    Straight synthesis with all-pass filter design based on
%	TEMPO analysis result
%   sy=straightSynthTB07ca(n2sgram,f0raw,f0var,f0varL,shiftm,fs, ...
%               pcnv,fconv,sconv,gdbw,delfrac,delsp,cornf,delfracind,ap,imap,imgi));
%	sy	: synthsized speech
%	n2sgram	: amplitude spectrogram
%	f0raw	: pitch pattern (Hz)
%	f0var	: expected F0 variation with fricative modification
%	f0varL	: expected F0 variation
%	shiftm	: frame shift (ms) for spectrogram
%	fs	: sampling freqnency (Hz)
%	pcnv	: pitch stretch factor
%	fconv	: freqnency stretch factor
%	sconv	: speaking duratin stretch factor (overridden if || imap || >1 )
%	gdbw	: finest resolution in group delay (Hz)
%	delfrac	: ratio of standard deviation of group delay in terms of F0
%	delsp	: standard deviation of group delay (ms)
%	cornf	: lower corner frequency for phase randomization (Hz)
%	delfracind	: selector of fixed and proportional group delay
%	ap	: aperiodicity measure
%	imap	: arbirtary mapping from new time (sample) to old time (frame)
%	imgi	: display indicator, 1: display on (default), 0: off
%   lowestF0    : lower limit of the resynthesized fundamental frequency (Hz)

%    Straight synthesis with all-pass filter design
%    by Hideki Kawahara
%	(c) ATR Human Info. Proc. Res. Labs. 1996
%	07/July/1996
%	12/Aug./1996
%	22/Aug./1996
%	06/Sep./1996 BUG FIX!!! wrong sign
%	07/Sep./1996 converted to function script
%	09/Sep./1996 coarse F0 information is possible
%	16/Sep./1996 tolerant to F0 extraction errors
%	02/Nov./1996 Now pitch extraction is perfect. No need for the hack.
%	02/Feb./1997 Without V/UV discrimination
%   08/June/1999 minor bug fix

%	03/Sep./1999 Graded excitation with one parameter
%	29/Nov./1999 Arbitrary time axis mapping
%	30/May/2001 revised aperiodicity control
%	08/April/2002 revised to remove magical LPF
%	11/August/2002 bug fix for V/UV transition
%	24/August/2002 more precise F0 control
% 	23/Sept./2002 minor adjustment for the length of the resynthesized signal
%	05/Dec./2002 minor bug fix based on M. Tsuzaki's comment
%   17/Dec./2002 bug fix in mid point selection
%	10/Aug./2005 modified by Takahashi on waitbar
%	10/Sept./2005 modified by Kawahara on waitbar
%   27/Nov./2005 modified by Kawahara for
%   21/April/2010 bug fix by Hideki Kawahara for aperiodicity
% 03/July/2016 refactored for MATLAB R2016a and Octave 4.0.2

%if nargin<=14; imgi=1; end; % 10/Sept./2005
statusReport = 'ok';% 27/Nov./2005
switch nargin % 27/Nov./2005
  case {1,2,3,4,5,6,7,8,9,10,11,12,13,14}
    imgi = 1;
    lowestF0 = 50;
  case {15}
    lowestF0 = 50;
end;
f0l=f0raw;
[nii,njj]=size(n2sgram);
njj=min([njj,length(f0raw)]);  % 18/Sep./1999
f0l=f0l(1:njj);  %03/Sep./1999
if min(f0l(f0l>0))*pcnv < lowestF0
  statusReport = ['Minimum synthesized F0 exceeded the lower limit(' num2str(lowestF0) ' Hz).'];
end;

fftLengthForLowestF0 = 2^ceil(log2(2*round(fs/lowestF0)));% 27/Nov./2005
fftl=nii+nii-2;
if fftl < fftLengthForLowestF0 % 27/Nov./2005
  niiNew = fftLengthForLowestF0/2+1;
  statusReport = 'The FFT length was inconsistent and replaced';
  n2sgram = interp1(0:nii-1,n2sgram,(0:niiNew-1)*(nii-1)/(niiNew-1));
  ap = interp1(0:nii-1,ap,(0:niiNew-1)*(nii-1)/(niiNew-1));
  fftl = fftLengthForLowestF0;
  nii = niiNew;
end;

% safeguard for ap mismatch 21/April/2010
if size(ap,1) ~= size(n2sgram,1)
  apDouble = zeros(size(n2sgram,1),size(ap,2));
  for ik = 1:size(ap,2)
    apDouble(:,ik) = interp1((0:size(ap,1)-1),ap(:,ik),...
      (0:size(n2sgram,1)-1)/((size(n2sgram,1)-1)/(size(ap,1)-1)),'linear','extrap');
  end;
  ap = apDouble;
end;

aprms=10.0.^(ap/20); % 23/Sept./1999
aprm=min(1,max(0.001,aprms*1.6-0.015)); % 30/May/2001

if length(fconv)==1
  idcv=min((0:fftl/2)/fconv+1,fftl/2+1); % f. stretch conv. tabel
elseif length(fconv)==nii
  idcv=fconv(:)';
elseif length(fconv) ~= nii
  idcv = 1:fftl/2+1;
  statusReport = [statusReport '\n' 'Frequency axix mapping function is not consistent with lowestF0.'];
end;
if length(imap)>1
  sy=zeros(length(imap)+3*fftl,1);disp('here!!');
else
  sy=zeros([round((njj*shiftm/1000*fs)*sconv+3*fftl+1),1]);
  imap=1:length(sy);
  imap=min(length(f0l),((imap-1)/fs*1000/shiftm/sconv+1));
end;
imap=[imap ones(1,round(fs*0.2))*length(f0l)]; % safe guard
ix=find(imap>=length(f0l), 1, 'first');
rmap=interp1(imap(1:ix),1:ix,1:length(f0l));

phs=fractpitch2(fftl); %  phs will have smooth phase function for unit delay

fftl2=fftl/2;
nsyn=length(sy);
idx=1;
bb=1:fftl;
rbb2=fftl/2:-1:2;

%------- parameters for noize based apf design
t=((1:fftl)-fftl/2-1)/fftl*2;
adjd=1.0./(1+exp(-20*t)); % correction function for smooth transition at fs/2
gw=exp(-0.25*pi*(fs*(t/2)/gdbw).^2); % slope difinition function
gw=gw/sum(gw);            % gdbw is the equvalent rectangular band width
fgw=real(fft(fftshift(gw))); % gw is the spectral smoothing window
df=fs/fftl*2*pi;    % normalization constant for integration and differentiation
fw=(1:fftl2+1)/fftl*fs; % frequency axis

trbw=300;		% width of transition area
rho=1.0./(1+exp(-(fw-cornf)/trbw)); % rondom group delay weighting function

%--------- frozen group delay component calculation ------
nz=randn(1,fftl2+1).*rho; % This is not effective. Left for randn status.
%---------
lft=1-hanning(fftl)+nz(1)*0; % +nz(1)*0 is dummy
lft=1.0./(1+exp(-(lft-0.5)*60));
ww=1.0./(1+exp(-(hanning(fftl)-0.3)*23)); % lifter for

iin=1;
if imgi==1; hpg=waitbar(0,'voiced part synthesis'); end; % 10/Aug./2005
icntr=0;
dmx=max(max(n2sgram));
while (idx < nsyn-fftl-10) && (ceil(iin)<length(f0l));
  
  icntr=icntr+1;
  iix=round(imap(round(idx)));
  ii=min(min(max(1,iix),njj),length(f0l));
  if imgi==1 && rem(icntr,10)==0;waitbar(ii/length(f0l));end; % 10/Aug./2005
  f0=f0l(round(ii));
  if f0==0
    f0=200;
  else
    f0=max(lowestF0/pcnv,f0l(round(ii))); % 27/Nov./2005
  end;
  f0=f0*pcnv;
  %- --------
  %	look ahead correction of F0 (cmmt by Simon Makin of Sheffield Univ.) 24/Augst/2002
  %- --------
  tnf0=fs/f0;
  tidx=idx+tnf0;
  tiix=round(imap(round(tidx)));
  tii=min(min(max(1,tiix),njj),length(f0l));
  tf0=f0l(round(tii));
  if (tf0>0) && (f0l(round(ii))>0)
    if f0l(round((ii+tii)/2))>0  % fix by H.K. on 17/Dec./2002
      f0=max(lowestF0/pcnv,f0l(round((ii+tii)/2))); % mid point
    else
      f0=f0l(round(ii));
    end;
    f0=f0*pcnv;
  end;
  
  %- --------
  ff=[n2sgram(round(idcv(:)),round(ii)); ...
    n2sgram(round(idcv(rbb2)),round(ii))];
  ccp=real(fft(log(ff+dmx/1000000)));  % 24 Sept. 1999 10000 -> 1000000
  ccp2=[ccp(1);2*ccp(2:fftl/2);0*ccp(fftl/2+1:fftl)];
  ffx=(fft(ccp2.*lft)/fftl);
  nidx=round(idx);
  
  nf0=fs/f0;
  frt=idx-nidx;
  frtz=exp(1i*phs*frt)'; % This was in a wrong sign!
  
  nz=randn(1,fftl2+1).*rho; %((rho*0+1)*mixNlow(round(ii))+(1-mixNlow(round(ii)))*rho);
  nz=real(ifft(fft([nz,nz(rbb2)]).*fgw));
  nz=nz*sqrt(fftl*gdbw/fs);   % correction factor for noise
  if delfracind, delsp=delfrac*1000/f0; end;
  nz=nz*delsp*df/1000;
  mz=cumsum([nz(1:fftl2+1),nz(rbb2)])-nz(1);
  mmz=-(mz-adjd*(rem((mz(fftl)+mz(2)),2*pi)-2*pi));
  pzr=exp(-1i*mmz)'; %.*[wlcut wlcut(rbb2)]'; % set ineffective 01/June/2001
  
  pz=pzr; % This makes random group delay to be effective
  wnz=aprm(round(idcv(:)),round(ii));  % 06/May/2001 This is correct!
  wpr=sqrt(max(0,1-wnz.*wnz));  % 23/Sept./1999
  
  rx=randn(round(nf0),1);
  %----------- temporal envelope control of the aperiodic component ---
  zt0=nf0/fs+rx(1)*0; % +rx(1)*0 is a dummy
  ztc=0.01; % time constant 10ms (for example)
  ztp=((1:round(nf0))'-1)/fs;
  nev=sqrt(2*zt0/ztc/(1-exp(-2*zt0/ztc)))*exp(-ztp/ztc);
  rx=randn(round(nf0),1);
  wfv=fft((rx-mean(rx)).*nev,fftl); % DC component removal 8/April/2002
  %--------------------------------------------------------------------
  
  ep=0*real(ffx);
  nf0n=round(nf0);
  gh=hanning(nf0n*2);
  ep(1:nf0n)=gh(nf0n:-1:1);
  ep(end:-1:end-nf0n+2)=ep(2:nf0n); % bug fix on 29/Jan./2003
  ep=-ep/sum(ep);
  ep(1)=ep(1)+1;
  epf=fft(ep);
  tx=fftshift(real(ifft(epf.*exp(ffx).*pz.*frtz.*[wpr;wpr(rbb2)]))).*ww; % 8/April/2002
  tx2=fftshift(real(ifft(exp(ffx).*frtz.*[wnz;wnz(rbb2)].*wfv))).*ww; % 31/May/2001
  sy(bb+nidx)=sy(bb+nidx)+(tx*sqrt(nf0)+tx2)*(f0raw(round(ii))>0); % 02/ Sept./1999
  idx=idx+nf0;
  iin=min(max(1,round(imap(round(idx)))),min(njj,length(f0raw))); % modification  on 5/Dec/2002 based on comments by M. Tsuzaki
  if (f0raw(round(ii))==0) && (f0raw(round(iin))>0) % (mixVlow(round(ii))<0.8) & (mixVlow(round(iin))>0.8)
    idxo=idx;
    ipos=find(f0raw(round(ii:iin))>0, 1, 'first')-1+ii;
    if isempty(ipos)
      idx=idxo;
    else
      idx=max(idxo-nf0+1,rmap(round(ipos))); % 11/August/2002 (Was -1 mistake??)
    end;
  end;
end;
if imgi==1; close(hpg); end; % 10/Aug./2005
ii=1;
idx=1;
f0=1000;

if imgi==1; hpg=waitbar(0,'unvoiced part synthesis'); end; % 10/Aug./2005
icntr=0;
while (idx < nsyn-fftl) && (ii<length(f0l));
  icntr=icntr+1;
  if imgi==1 && rem(icntr,10)==0;waitbar(ii/length(f0l));end;
  nidx=round(idx);
  if f0raw(ii) == 0; %mixNhigh(ii)>0.03
    ff=[n2sgram(round(idcv(:)),ii);n2sgram(round(idcv(rbb2)),ii)];
    ccp=real(fft(log(ff+dmx/100000))); % 23rd July, 1999 % 24th Sept.
    ccp2=[ccp(1);2*ccp(2:fftl/2);0*ccp(fftl/2+1:fftl)];
    ffx=(fft(ccp2.*lft)/fftl);
    nf0=fs/f0;
    tx=fftshift(real(ifft(exp(ffx))));
    rx=randn([round(nf0),1]);
    tnx=fftfilt(rx-mean(rx),tx); % DC component removal 8/April/2002
    sy(bb+nidx)=sy(bb+nidx)+tnx(bb).*ww;
  end;
  idx=idx+nf0;
  ii=round(imap(round(idx)));
end;
if imgi==1; close(hpg); end; % 10/Aug./2005
sy2=sy(fftl/2+(1:ix));
sy=sy2;
switch nargout
  case {1}
  case {2}
    synthSataus = statusReport;
end;
end
