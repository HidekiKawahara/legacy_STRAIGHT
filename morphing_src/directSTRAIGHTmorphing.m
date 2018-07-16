function mObject3 = directSTRAIGHTmorphing(mObject1,mObject2,mRate,mixMethod);
%   Morphing based on direct mixing of STRAIGHT parameters
%   (Without time alignment)
%   mObject3 = directSTRAIGHTmorphing(mObject1,mObject2,mRate);

%   Designed and coded by Hideki Kawahara
%   27/Feb./2005
%   Copyright(c) 2005, Hideki Kawahara

if mObject1.samplingFrequency ~= mObject2.samplingFrequency
    mObject3 = [];
    return
end;
if mObject1.frameUpdateInterval ~= mObject2.frameUpdateInterval
    mObject3 = [];
    return
end;

nw1 = length(mObject1.F0);
nw2 = length(mObject2.F0);
[nr1,nc1] = size(mObject1.spectrogram);
[nr2,nc2] = size(mObject2.spectrogram);
nr3 = max(nr1,nr2);
nc3 = max(max(nc1,nc2),max(nw1,nw2));
nsg = zeros(nr3,nc3);
nSgram = zeros(nr3,nc3);
ap = zeros(nr3,nc3);
f0 = zeros(nc3);
nVoice = zeros(nc3);

switch mixMethod
    case 'linear'
        nsg(1:nr1,1:nc1) = (1-mRate)*mObject1.spectrogram;
        nSgram(1:nr1,1:nc1) = nSgram(1:nr1,1:nc1)+(1-mRate);
        nsg(1:nr2,1:nc2) = mRate*mObject2.spectrogram+nsg(1:nr2,1:nc2);
        nSgram(1:nr2,1:nc2) = nSgram(1:nr2,1:nc2)+mRate;
        nsg = nsg./nSgram;
    case 'log'
        nsg(1:nr1,1:nc1) = (1-mRate)*log(mObject1.spectrogram);
        nSgram(1:nr1,1:nc1) = nSgram(1:nr1,1:nc1)+(1-mRate);
        nsg(1:nr2,1:nc2) = mRate*log(mObject2.spectrogram)+nsg(1:nr2,1:nc2);
        nSgram(1:nr2,1:nc2) = nSgram(1:nr2,1:nc2)+mRate;
        nsg = exp(nsg./nSgram);
end;
ap(1:nr1,1:nc1) = (1-mRate)*mObject1.aperiodicityIndex;
ap(1:nr2,1:nc2) = mRate*mObject2.aperiodicityIndex+ap(1:nr2,1:nc2);

f0(mObject1.F0>0) = (1-mRate)*log(mObject1.F0(mObject1.F0>0));
nVoice(mObject1.F0>0) = nVoice(mObject1.F0>0)+(1-mRate);
f0(mObject2.F0>0) = mRate*log(mObject2.F0(mObject2.F0>0))+f0(mObject2.F0>0);
nVoice(mObject2.F0>0) = nVoice(mObject2.F0>0)+mRate;
f0(nVoice>0) = exp(f0(nVoice>0)./nVoice(nVoice>0));

mObject3 = createMobject;
mObject3 = updateFieldOfMobject(mObject3,'spectrogram',nsg);
mObject3 = updateFieldOfMobject(mObject3,'aperiodicityIndex',ap);
mObject3 = updateFieldOfMobject(mObject3,'F0',f0);
