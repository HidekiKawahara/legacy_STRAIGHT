function sy=straightSynthTB06(n2sgram,f0raw,f0var,f0varL,shiftm,fs, ...
            pcnv,fconv,sconv,gdbw,delfrac,delsp,cornf,delfracind);
%    Straight synthesis with all-pass filter design based on
%	TEMPO analysis result
%   sy=straightSynthTB06(n2sgram,f0raw,f0var,f0varL,shiftm,fs, ...
%               pcnv,fconv,sconv,gdbw,delfrac,delsp,cornf,delfracind);
%	sy	: synthsized speech
%	n2sgram	: amplitude spectrogram
%	f0raw	: pitch pattern (Hz) 
%	f0var	: expected F0 variation with fricative modification
%	f0varL	: expected F0 variation
%	shiftm	: frame shift (ms) for spectrogram
%	fs	: sampling freqnency (Hz)
%	pcnv	: pitch stretch factor
%	fconv	: freqnency stretch factor 
%	sconv	: speaking duratin stretch factor
%	gdbw	: finest resolution in group delay (Hz)
%	delfrac	: ratio of standard deviation of group delay in terms of F0
%	delsp	: standard deviation of group delay (ms)
%	cornf	: lower corner frequency for phase randomization (Hz)
%	delfracind	: selector of fixed and proportional group delay


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

f0l=f0raw;
[nii,njj]=size(n2sgram);
fftl=nii+nii-2;
fftl2=fftl/2;
if length(fconv)==1
    idcv=min([0:fftl/2]/fconv+1,fftl/2+1); % f. stretch conv. tabel
  elseif length(fconv)==nii
	idcv=fconv(:)';
end;
sy=zeros([round((njj*shiftm/1000*fs)*sconv+3*fftl+1),1]);
syo=sy;

mixVhigh=sqrt(0.25./(f0var+0.25));
mixNhigh=sqrt(1-0.25./(f0var+0.25));

mixVlow=sqrt(0.25./(f0varL+0.25));
mixNlow=sqrt(1-0.25./(f0varL+0.25));

phs=fractpitch2(fftl); %  phs will have smooth phase function for unit delay

a=([0:fftl2-1,0,-(fftl2-1:-1:1)])/fftl2;
sz=a'*pi;
ta=[0:fftl2-1]/fftl2/2*2*pi;
t=[ta,0,-ta(fftl2:-1:2)];

fftl2=fftl/2;
nsyn=length(sy);
idx=1;
bb=1:fftl;
bb2=1:fftl2;
rbb2=fftl/2:-1:2;

%------- shaping for low-frequency noize supression

fxa=(0:fftl2)/fftl*fs;
f0tmp=f0l.*(mixVlow>0.8);
lowcutf=mean(f0tmp(f0tmp>0))*0.7*pcnv;
%lowcutfav=mean(f0l(f0l>0))*0.8;
lowcutfav=lowcutf;
%wlcutav=1.0./(1+exp(-5*(fxa-lowcutfav)/(lowcutfav/3)));
wlcutav=1.0./(1+exp(-14*(fxa-lowcutfav)/(lowcutfav/1)));
%keyboard;


%------- parameters for noize based apf design

t=([1:fftl]-fftl/2-1)/fftl*2;
adjd=1.0./(1+exp(-20*t)); % correction function for smooth transition at fs/2
gw=exp(-0.25*pi*(fs*(t/2)/gdbw).^2); % slope difinition function
gw=gw/sum(gw);            % gdbw is the equvalent rectangular band width
fgw=real(fft(fftshift(gw))); % gw is the spectral smoothing window
df=fs/fftl*2*pi;    % normalization constant for integration and differentiation
fw=(1:fftl2+1)/fftl*fs; % frequency axis

trbw=300;		% width of transition area
rho=1.0./(1+exp(-(fw-cornf)/trbw)); % rondom group delay weighting function

[snn,smm]=size(n2sgram);
fqx=(0:snn-1)/snn*fs/2;
chigh=1.0./(1+exp(-(fqx-600)/100))';
clow=1.0-chigh;

f0arc=0;
lft=1-hanning(fftl);
lft=1.0./(1+exp(-(lft-0.5)*60));
ww=1.0./(1+exp(-(hanning(fftl)-0.3)*23)); % lifter for

iin=1;
dmx=max(max(n2sgram));
while (idx < nsyn-fftl-10) & (ceil(iin)<length(f0l));
  
  iix=idx/fs*1000/shiftm/sconv+1;
  ii=min(min(max(1,iix),njj),length(f0l));
%	wlcut=wlcutav;
    f0=max(40,f0l(round(ii)));
    f0=f0*pcnv;
%    wlcut=1.0./(1+exp(-7*(fxa-f0*0.9)/lowcutf));  % 23rd July, 1999
    wlcut=1.0./(1+exp(-10*(fxa-f0*0.7)/f0));  % 23rd July, 1999

  mix=mixVlow(round(ii))*clow(round(idcv(:)))+mixVhigh(round(ii))*chigh(round(idcv(:)));
  ff=[n2sgram(round(idcv(:)),round(ii)); ...
      n2sgram(round(idcv(rbb2)),round(ii))];
  ff=ff.*[wlcut wlcut(rbb2)]';
  ccp=real(fft(log(ff+dmx/1000000))); % 24 Sept. 1999  10000 -> 1000000
  ccp2=[ccp(1);2*ccp(2:fftl/2);0*ccp(fftl/2+1:fftl)];
  ffx=(fft(ccp2.*lft)/fftl);
  nidx=round(idx);


%    wlcut=1.0./(1+exp(-20*(fxa-lowcutf)/lowcutf));
    nf0=fs/f0;
    frt=idx-nidx;
    frtz=exp(i*phs*frt)'; % This was in a wrong sign!
    
    nz=randn(1,fftl2+1).*((rho*0+1)*mixNlow(round(ii))+(1-mixNlow(round(ii)))*rho);
    nz=real(ifft(fft([nz,nz(rbb2)]).*fgw));
    nz=nz*sqrt(fftl*gdbw/fs);   % correction factor for noise
    if delfracind, delsp=delfrac*1000/f0; end;
    nz=nz*delsp*df/1000;
    mz=cumsum([nz(1:fftl2+1),nz(rbb2)])-nz(1);
    mmz=-(mz-adjd*(rem((mz(fftl)+mz(2)),2*pi)-2*pi));
    pz=exp(-i*mmz)'; %.*[wlcut wlcut(rbb2)]';

    tx=fftshift(real(ifft(exp(ffx).*pz.*frtz.*[mix;mix(rbb2)]))).*ww;
%    tx=fftshift(real(ifft(ff.*pz.*frtz.*[mix;mix(rbb2)]))).*ww;
    sy(bb+nidx)=sy(bb+nidx)+tx*sqrt(nf0);
%	if abs(round(ii)-90)<10
%		keyboard;
%	end;	

  idx=idx+nf0;
  iin=min(length(f0l),idx/fs*1000/shiftm/sconv+1);
  if (mixVlow(round(ii))<0.8) & (mixVlow(round(iin))>0.8)
    idxo=idx;
    ipos=min(find(mixVlow(round(ii:iin))>0.8))-1+ii;
    if length(ipos)==0
      idx=idxo;
    else
      idx=max(idxo-nf0+1,(ipos-1)*fs/1000*shiftm*sconv);
    end;
  end;
%  disp([idx,iin])
end;
%sy=sy*0;
ii=1;
idx=1;
f0=500;
f0=1000;
%wlcutfric=1.0./(1+exp(-14*(fxa-lowcutfav*2)/(lowcutfav)));
wlcutfric=1.0./(1+exp(-14*(fxa-lowcutfav)/(lowcutfav))); % 31/July/1999

while (idx < nsyn-fftl) & (ii<length(f0l));
  ii=round(min([length(f0l),idx/fs*1000/shiftm/sconv+1,njj]));
  nidx=round(idx);
  if mixNhigh(ii)>0.03
    mix=mixNlow(ii)*clow(round(idcv(:)))+mixNhigh(ii)*chigh(round(idcv(:)));
    ff=[n2sgram(round(idcv(:)),ii);n2sgram(round(idcv(rbb2)),ii)];
%    ff=ff.*[wlcut wlcut(rbb2)]';
    ff=ff.*[wlcutfric wlcutfric(rbb2)]';
%    ccp=real(fft(log(ff+0.001))); % 23rd July, 1999
	ccp=real(fft(log(ff+dmx/100000))); % 23rd July, 1999 % 24th Sept. 1999
    ccp2=[ccp(1);2*ccp(2:fftl/2);0*ccp(fftl/2+1:fftl)];
    ffx=(fft(ccp2.*lft)/fftl);
    nf0=fs/f0;

%============= deleted on 18/July/1999 ======
%    if f0l(ii) > 0
%      f0x=lowcutf; % f0l(ii)*pcnv;
%      f0x=f0l(ii)*pcnv;
%      wlcut=1.0./(1+exp(-20*(fxa-f0x*0.8)/lowcutf));
%	  wlcut=wlcutav;
%      tx=fftshift(real(ifft(exp(ffx).*[wlcut.*mix' wlcut(rbb2).*mix(rbb2)']')));
%    else
%      tx=fftshift(real(ifft(exp(ffx).*[wlcutav.*mix' wlcutav(rbb2).*mix(rbb2)']')));
%    end;
	tx=fftshift(real(ifft(exp(ffx))));
%============= end of modification on 18/July/1999 ====

    rx=randn([round(nf0),1]);
    tnx=fftfilt(rx,tx);
    sy(bb+nidx)=sy(bb+nidx)+tnx(bb).*ww;
  end;
  idx=idx+nf0;
  ii=min(length(f0l),idx/fs*1000/shiftm/sconv+1);
end;

sy2=sy(fftl/2+(1:round((njj*shiftm/1000*fs)*sconv)));

lowcutf=70;
if lowcutf <70
  lowcutf=70;
end;
%[b,a]=butter(5,lowcutf/fs*2,'high');
%sy=filter(b,a,sy2);
sy=sy2;


