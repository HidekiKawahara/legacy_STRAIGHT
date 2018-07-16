function mObject3 = timeAlignedDirectSTRAIGHTmorphing(mObject1,mObject2,mRate,mixMethod);
%   Morphing based on time-aligned mixing of STRAIGHT parameters
%   mObject3 = timeAlignedDirectSTRAIGHTmorphing(mObject1,mObject2,mRate,mixMethod);

%   Designed and coded by Hideki Kawahara
%   28/Feb./2005
%   Copyright(c) 2005, Hideki Kawahara

mObject3 = checkForSimilarity(mObject1,mObject2);
if length(mObject3) ==0;return;end;
dtFrame = mObject1.frameUpdateInterval;
endLocation1 = (length(mObject1.F0)-1)*dtFrame; % in ms
endLocation2 = (length(mObject2.F0)-1)*dtFrame; % in ms
timeAnchor1 = [0;mObject1.anchorTimeLocation;endLocation1];
timeAnchor2 = [0;mObject2.anchorTimeLocation;endLocation2];
locationOn1 = (0:length(mObject1.F0)-1)*dtFrame;
locationOn2 = (0:length(mObject2.F0)-1)*dtFrame;
mapFrom1to2 = interp1(timeAnchor1,timeAnchor2,locationOn1);
[nr1,nc1] = size(mObject1.spectrogram);
[nr2,nc2] = size(mObject2.spectrogram);

%---- mixing on mObject1's time axis
nAxis1 = length(locationOn1); 
nAxis2 = length(locationOn2);
morphedF0 = zeros(nAxis1,1);
morphedAp = zeros(nr1,nAxis1);
morphedSgram = zeros(nr1,nAxis1);
weightSumF0 = zeros(nAxis1,1);
for ii=1:nAxis1
    mappedIndexOn2 = mapFrom1to2(ii)/dtFrame+1;
    iFloor = floor(mappedIndexOn2);
    iFraction = mappedIndexOn2-iFloor;
    dAp = iFraction*(mObject2.aperiodicityIndex(:,min(iFloor+1,nAxis2))-mObject2.aperiodicityIndex(:,iFloor));
    morphedAp(:,ii) = (1-mRate)*mObject1.aperiodicityIndex(:,ii)+mRate*(mObject2.aperiodicityIndex(:,iFloor)+dAp);
    switch mixMethod
        case 'linear'
            dSgram = iFraction*(mObject2.spectrogram(:,min(iFloor+1,nAxis2))-mObject2.spectrogram(:,iFloor));
            morphedSgram(:,ii) = (1-mRate)*mObject1.spectrogram(:,ii)+mRate*(mObject2.spectrogram(:,iFloor)+dSgram);
        case 'log'
            dSgram = iFraction*(log(mObject2.spectrogram(:,min(iFloor+1,nAxis2)))-log(mObject2.spectrogram(:,iFloor)));
            tmp = (1-mRate)*log(mObject1.spectrogram(:,ii))+mRate*(log(mObject2.spectrogram(:,iFloor))+dSgram);
            morphedSgram(:,ii) = exp(tmp);
    end;
    if mObject1.F0(ii)>0
        morphedF0(ii) = (1-mRate)*log(mObject1.F0(ii));
        weightSumF0(ii) = (1-mRate);
    end;
    if (mObject2.F0(iFloor)>0) & (mObject2.F0(min(iFloor+1,nAxis2))>0)
        dF0 = iFraction*(log(mObject2.F0(min(iFloor+1,nAxis2)))-log(mObject2.F0(iFloor)));
        morphedF0(ii) = mRate*(log(mObject2.F0(iFloor))+dF0)+morphedF0(ii);
        weightSumF0(ii) = weightSumF0(ii)+mRate;
    end;
end;
morphedF0(weightSumF0>0) = exp(morphedF0(weightSumF0>0)./weightSumF0(weightSumF0>0));

%----- mapping back onto morphed time axis
timeAnchorMorph = (1-mRate)*timeAnchor1 + mRate*timeAnchor2;
locationOnMorph = (0:(timeAnchorMorph(end)/dtFrame))*dtFrame;
mapFormMorphTo1 = interp1(timeAnchorMorph,timeAnchor1,locationOnMorph);
nAxisMorph = length(locationOnMorph);
morphedApOnMorph = zeros(nr1,nAxisMorph);
morphedSgramOnMorph = zeros(nr1,nAxisMorph);
morphedF0onMorph = zeros(nAxisMorph,1);
for ii=1:nAxisMorph
    mappedIndexOn1 = mapFormMorphTo1(ii)/dtFrame+1;
    iFloor = floor(mappedIndexOn1);
    iFraction = mappedIndexOn1-iFloor;
    morphedApOnMorph(:,ii) = morphedAp(:,iFloor) ...
        +iFraction*(morphedAp(:,min(iFloor+1,nAxis1))-morphedAp(:,iFloor));
    morphedSgramOnMorph(:,ii) = morphedSgram(:,iFloor) ...
        +iFraction*(morphedSgram(:,min(iFloor+1,nAxis1))-morphedSgram(:,iFloor));
    if (morphedF0(iFloor)>0) & (morphedF0(min(iFloor+1,nAxis1))>0)
        dF0 = iFraction*(morphedF0(min(iFloor+1,nAxis1))-morphedF0(iFloor));
        morphedF0onMorph(ii) = morphedF0(iFloor)+dF0;
    end;
end;
mObject3.F0 = morphedF0onMorph; 
mObject3.aperiodicityIndex = morphedApOnMorph;
mObject3.spectrogram = morphedSgramOnMorph;
mObject3.anchorTimeLocation = timeAnchorMorph(2:end-1);
mObject3.anchorFrequency = (1-mRate)*mObject1.anchorFrequency+mRate*mObject2.anchorFrequency;
%mObject3 = morphedAp; % This line is a dummy.

%%% ------ Internal function to check for object's similarity
function mObject3 = checkForSimilarity(mObject1,mObject2)
mObject3 = [];
if mObject1.samplingFrequency ~= mObject2.samplingFrequency;mObject3 = [];return;end;
if mObject1.frameUpdateInterval ~= mObject2.frameUpdateInterval;mObject3 = [];return;end;
if length(mObject1.anchorTimeLocation) ~= length(mObject2.anchorTimeLocation);mObject3 = [];return;end;
nAnchor = length(mObject1.anchorTimeLocation);
for ii=1:nAnchor % check for similarity of anchor structure
    frequencyAnchor1 = mObject1.anchorFrequency;
    frequencyAnchor2 = mObject2.anchorFrequency;
    if (sum(frequencyAnchor1>0) ~= sum(frequencyAnchor2>0)) | ...
            (sum(frequencyAnchor1<0) ~= sum(frequencyAnchor2<0))
        return;
    end;
end;
mObject3 = createMobject;
