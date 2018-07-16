function [ap,analysisParams]=exstraightAPind(x,fs,f0,optionalParams)
%   Aperiodicity index extraction for STRAIGHT
%   [ap,analysisParams]=exstraightAPind(x,fs,f0,optionalParams)
%   Input parameters
%   x   : input signal. if it is multi channel, only the first channel is used
%   fs  : sampling frequency (Hz)
%   f0  : fundamental frequency (Hz)
%   optionalParams : Optional parameters for analysis
%   Output parameters
%   ap  : amount of aperiodic component in the time frequency represntation
%       : represented in dB
%   analysisParams : Analysis parameters actually used
%
%   Usage:
%   Case 1: The simplest method
%   ap=exstraightAPind(x,fs,f0); 
%   Case 2: You can get to know what parameters were used.
%   [ap,analysisParams]=exstraightAPind(x,fs,f0);
%   CAse 3: You can have full control of STRAIGHT synthesis.
%       Please use case 2 to find desired parameters to modify.
%   [ap,analysisParams]=exstraightAPind(x,fs,f0,optionalParams);

%   Notes on programing style
%   This routine is based on the current (2005.1.31) implementation of
%   STRAIGHT that consist of many legacy fragments. They were intentionally
%   kept for maintaining historic record. Revised functions written in a
%   reasonable stylistic practice will be made available soon.

%   Designed and coded by Hideki Kawahara
%   15/January/2005
%   01/February/2005 extended for user control
%   13/March/2005 Aperiodicity index extraction part is isolated
%	30/April/2005 modification for Matlab v7.0 compatibility
%   11/Sept./2005 waitbar control is fixed.
%   05/July/2006 default values are modified, framem
		
%---Check for number of input parameters
switch nargin
    case 3
        prm=zinitializeParameters;
    case 4
        prm=replaceSuppliedParameters(optionalParams);
    otherwise
        disp('Number of arguments is 3 or 4!');
        return;
end

%   Initialize default parameters
f0ceil = prm.F0searchUpperBound; % f0ceil
framem = prm.F0defaultWindowLength; % default frame length for pitch extraction (ms)
f0shiftm = prm.F0frameUpdateInterval; % shiftm % F0 calculation interval (ms)

fftl=1024;	% default FFT length

framel=framem*fs/1000;

if fftl < framel
    fftl=2^ceil(log(framel)/log(2));
end;

[nr,nc]=size(x);
if nr>nc
    x=x(:,1);
else
    x=x(1,:)';
end;

imageOn = prm.DisplayPlots; % imgi=1; % image display indicator (1: display image)

% paramaters for F0 refinement
fftlf0r = prm.refineFftLength; %fftlf0r=1024; % FFT length for F0 refinement
tstretch = prm.refineTimeStretchingFactor; %tstretch=1.1; % time window stretching factor
nhmx = prm.refineNumberofHarmonicComponent; %nhmx=3; % number of harmonic components for F0 refinement
iPeriodicityInterval = prm.periodicityFrameUpdateInterval; % frame update interval for periodicity index (ms)

%---- F0 refinement 
nstp=1; % start position of F0 refinement (samples)
nedp=length(f0); % last position of F0 refinement (samples)
dn=floor(fs/(f0ceil*3*2)); % fix by H.K. at 28/Jan./2003
[f0raw,ecr]=refineF06(decimate(x,dn),fs/dn,f0,fftlf0r,tstretch,nhmx,f0shiftm,nstp,nedp,imageOn); % 31/Aug./2004

ecrt=ecr;
ecrt(f0raw==0)=ecrt(f0raw==0)*NaN;

%----- aperiodicity estimation
f0raw=f0;
[apvq,dpvq,~,~]=aperiodicpartERB2(x,fs,f0raw,f0shiftm,iPeriodicityInterval,fftl/2+1,imageOn); % 10/April/2002
apv=10*log10(apvq); % for compatibility
dpv=10*log10(dpvq); % for compatibility

%- ---------
%   Notes on aperiodicity estimation: The previous implementation of
%   aperiodicity estimation was sensitive to low frequency noise. It is a
%   bad news, because environmental noise usually has its power in the low
%   frequency region. The following corrction uses the C/N information
%   which is the byproduct of fixed point based F0 estimation.
%   by H.K. 04/Feb./2003
%- ---------
dpv=correctdpv(apv,dpv,iPeriodicityInterval,f0raw,ecrt,f0shiftm,fs); % Aperiodicity correction 04/Feb./2003 by H.K.

if imageOn
    bv=boundmes2(apv,dpv,fs,f0shiftm,iPeriodicityInterval,fftl/2+1);
    figure;
    semilogy((0:length(bv)-1)*f0shiftm,0.5./10.0.^(bv));grid on;
    set(gcf,'PaperPosition', [0.634517 0.634517 19.715 28.4084]);
end;

ap=aperiodiccomp(apv,dpv,iPeriodicityInterval,f0raw,f0shiftm,imageOn); % 11/Sept./2005

switch nargout
    case 1
    case 2
        analysisParams=prm;
    otherwise
        disp('Number of output parameters has to be 1 or 2!')
end;
end

%%%---- internal functions

%%%------
function prm=zinitializeParameters
prm.F0searchLowerBound=40; % f0floor
prm.F0searchUpperBound=800; % f0ceil
prm.F0defaultWindowLength = 80; % default frame length for pitch extraction (ms)
prm.F0frameUpdateInterval=1; % shiftm % F0 calculation interval (ms)
prm.NofChannelsInOctave=24; % nvo=24; % Number of channels in one octave
prm.IFWindowStretch=1.2; % mu=1.2; % window stretch from isometric window
prm.DisplayPlots=0; % imgi=1; % image display indicator (1: display image)
prm.IFsmoothingLengthRelToFc=1; %  smp=1; % smoothing length relative to fc (ratio)
prm.IFminimumSmoothingLength=5; %  minm=5; % minimum smoothing length (ms)
prm.IFexponentForNonlinearSum=0.5; % pc=0.5; % exponent to represent nonlinear summation
prm.IFnumberOfHarmonicForInitialEstimate=1; % nc=1; % number of harmonic component to use (1,2,3)
prm.refineFftLength=1024; %fftlf0r=1024; % FFT length for F0 refinement
prm.refineTimeStretchingFactor=1.1; %tstretch=1.1; % time window stretching factor
prm.refineNumberofHarmonicComponent=3; %nhmx=3; % number of harmonic components for F0 refinement
prm.periodicityFrameUpdateInterval=5; % frame update interval for periodicity index (ms)return
prm.note=' '; % Any text to be printed on the source information plot
end

%%%--------
function prm=replaceSuppliedParameters(prmin)
prm=zinitializeParameters;
if isfield(prmin,'F0searchLowerBound')==1;
    prm.F0searchLowerBound=prmin.F0searchLowerBound;end;
if isfield(prmin,'F0searchUpperBound')==1;
    prm.F0searchUpperBound=prmin.F0searchUpperBound;end;
if isfield(prmin,'F0defaultWindowLength')==1;
    prm.F0defaultWindowLength=prmin.F0defaultWindowLength;end;
if isfield(prmin,'F0frameUpdateInterval')==1;
    prm.F0frameUpdateInterval=prmin.F0frameUpdateInterval;end;
if isfield(prmin,'NofChannelsInOctave')==1;
    prm.NofChannelsInOctave=prmin.NofChannelsInOctave;end;
if isfield(prmin,'IFWindowStretch')==1;
    prm.IFWindowStretch=prmin.IFWindowStretch;end;
if isfield(prmin,'DisplayPlots')==1;
    prm.DisplayPlots=prmin.DisplayPlots;end;
if isfield(prmin,'IFsmoothingLengthRelToFc')==1;
    prm.IFsmoothingLengthRelToFc=prmin.IFsmoothingLengthRelToFc;end;
if isfield(prmin,'IFminimumSmoothingLength')==1;
    prm.IFminimumSmoothingLength=prmin.IFminimumSmoothingLength;end;
if isfield(prmin,'IFexponentForNonlinearSum')==1;
    prm.IFexponentForNonlinearSum=prmin.IFexponentForNonlinearSum;end;
if isfield(prmin,'IFnumberOfHarmonicForInitialEstimate')==1;
    prm.IFnumberOfHarmonicForInitialEstimate=prmin.IFnumberOfHarmonicForInitialEstimate;end;
if isfield(prmin,'refineFftLength')==1;
    prm.refineFftLength=prmin.refineFftLength;end;
if isfield(prmin,'refineTimeStretchingFactor')==1;
    prm.refineTimeStretchingFactor=prmin.refineTimeStretchingFactor;end;
if isfield(prmin,'refineNumberofHarmonicComponent')==1;
    prm.refineNumberofHarmonicComponent=prmin.refineNumberofHarmonicComponent;end;
if isfield(prmin,'periodicityFrameUpdateInterval')==1;
    prm.periodicityFrameUpdateInterval=prmin.periodicityFrameUpdateInterval;end;
if isfield(prmin,'note')==1;
    prm.note=prmin.note;end;
end
