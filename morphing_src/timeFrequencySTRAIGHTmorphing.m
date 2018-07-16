function mObject3 = timeFrequencySTRAIGHTmorphing(mObject1,mObject2,mRate,mixMethod);
%   Morphing based on STRAIGHT parameters
%   mObject3 = timeFrequencySTRAIGHTmorphing(mObject1,mObject2,mRate,mixMethod);

%   Designed and coded by Hideki Kawahara
%   28/Feb./2005
%   Copyright(c) 2005, Hideki Kawahara
%   14/March/2005 bug fix on sampling frequency
%   01/Oct./2005 bug fix on similarity check
%   04/Oct./2005 partial morphing extension
%   18/Oct./2005 direct differential manipulation and API cange
%   29/Jan./2006 bug fix on boundary conditions

switch nargin
    case 0
        mObject3.morphingObject = createMobject;
        mixRate.F0 = 0;
        mixRate.spectrum = 0;
        mixRate.aperiodicity = 0;
        mixRate.coordinate = 0;
        mObject3.mixRate = mixRate;
        mObject3.mixMethods = {'linear','log','differentialLogarithm'};
        return
end;
mObject3 = checkForSimilarity(mObject1,mObject2);
mixRate = checkForMorphingConditions(mRate);
mObject1 = checkForIntegrity(mObject1);
mObject2 = checkForIntegrity(mObject2);
fs = mObject1.samplingFrequency;
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

%---- initialize frequency mapping function
fmapFrom1to2OnTime1 = generateFrequencyMap(mObject1,mObject2);

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
    fIndex = floor(fmapFrom1to2OnTime1(:,ii)/fs*2*(nr1-1))+1;
    dAp = iFraction*(mObject2.aperiodicityIndex(:,min(iFloor+1,nAxis2))-mObject2.aperiodicityIndex(:,min(iFloor,nAxis2)));
    ap2on2faxis = mObject2.aperiodicityIndex(:,min(iFloor,nAxis2))+dAp;
    ap2on1faxis = ap2on2faxis(fIndex);
    morphedAp(:,ii) = (1-mixRate.aperiodicity)*mObject1.aperiodicityIndex(:,ii)+mixRate.aperiodicity*ap2on1faxis; %04/Oct/2005 HK
    switch mixMethod
        case 'linear'
            dSgram = iFraction*(mObject2.spectrogram(:,min(iFloor+1,nAxis2))-mObject2.spectrogram(:,min(iFloor,nAxis2)));
            sgram2on2faxis = mObject2.spectrogram(:,min(iFloor,nAxis2))+dSgram;
            sgram2on1faxis = sgram2on2faxis(fIndex);
            morphedSgram(:,ii) = (1-mixRate.spectrum)*mObject1.spectrogram(:,ii)+mixRate.spectrum*sgram2on1faxis;
        case 'log'
            dSgram = iFraction*(log(mObject2.spectrogram(:,min(iFloor+1,nAxis2)))-log(mObject2.spectrogram(:,min(iFloor,nAxis2))));
            sgram2on2faxis = log(mObject2.spectrogram(:,min(iFloor,nAxis2)))+dSgram;
            sgram2on1faxis = sgram2on2faxis(fIndex);
            tmp = (1-mixRate.spectrum)*log(mObject1.spectrogram(:,ii))+mixRate.spectrum*sgram2on1faxis;
            morphedSgram(:,ii) = exp(tmp);
        case 'differentialLogarithm'
            dSgram = iFraction*(mObject2.spectrogram(:,min(iFloor+1,nAxis2))-mObject2.spectrogram(:,min(iFloor,nAxis2)));
            sgram2on2faxis = mObject2.spectrogram(:,min(iFloor,nAxis2))+dSgram;
            sgram2on1faxis = sgram2on2faxis(fIndex);
            tmp = (1-mixRate.spectrum)*log(mObject1.spectrogram(:,ii))+mixRate.spectrum*sgram2on1faxis;
            morphedSgram(:,ii) = exp(tmp);
    end;
    if mObject1.F0(ii)>0
        morphedF0(ii) = (1-mixRate.F0)*log(mObject1.F0(ii));
        weightSumF0(ii) = (1-mixRate.F0);
    end;
    if (mObject2.F0(iFloor)>0) & (mObject2.F0(min(iFloor+1,nAxis2))>0)
        dF0 = iFraction*(log(mObject2.F0(min(iFloor+1,nAxis2)))-log(mObject2.F0(min(iFloor,nAxis2))));
        morphedF0(ii) = mixRate.F0*(log(mObject2.F0(min(iFloor,nAxis2)))+dF0)+morphedF0(ii);
        weightSumF0(ii) = weightSumF0(ii)+mixRate.F0;
    end;
end;
morphedF0(weightSumF0>0) = exp(morphedF0(weightSumF0>0)./weightSumF0(weightSumF0>0));

%----- mapping back onto morphed time axis
timeAnchorMorph = (1-mixRate.coordinate)*timeAnchor1 + mixRate.coordinate*timeAnchor2;
locationOnMorph = (0:(timeAnchorMorph(end)/dtFrame))*dtFrame;
mapFormMorphTo1 = interp1(timeAnchorMorph,timeAnchor1,locationOnMorph);
nAxisMorph = length(locationOnMorph);
morphedApOnMorph = zeros(nr1,nAxisMorph);
morphedSgramOnMorph = zeros(nr1,nAxisMorph);
morphedF0onMorph = zeros(nAxisMorph,1);
%----- set place holders
mObject3.samplingFrequency = fs;
mObject3.F0 = morphedF0onMorph; 
mObject3.aperiodicityIndex = morphedApOnMorph;
mObject3.spectrogram = morphedSgramOnMorph;
mObject3.anchorTimeLocation = timeAnchorMorph(2:end-1);
mObject3.anchorFrequency = (1-mixRate.coordinate)*mObject1.anchorFrequency+mixRate.coordinate*mObject2.anchorFrequency;
%------ nitialize frequency mapping function
fmapFromMorphto1OnTimeMorph = generateFrequencyMap(mObject3,mObject1);
for ii=1:nAxisMorph
    mappedIndexOn1 = mapFormMorphTo1(ii)/dtFrame+1;
    iFloor = floor(mappedIndexOn1);
    iFraction = mappedIndexOn1-iFloor;
    fIndex = floor(fmapFromMorphto1OnTimeMorph(:,ii)/fs*2*(nr1-1))+1;
    morphedApOnMorph(:,ii) = morphedAp(fIndex,iFloor) ...
        +iFraction*(morphedAp(fIndex,min(iFloor+1,nAxis1))-morphedAp(fIndex,iFloor));
    morphedSgramOnMorph(:,ii) = morphedSgram(fIndex,iFloor) ...
        +iFraction*(morphedSgram(fIndex,min(iFloor+1,nAxis1))-morphedSgram(fIndex,iFloor));
    if (morphedF0(iFloor)>0) & (morphedF0(min(iFloor+1,nAxis1))>0)
        dF0 = iFraction*(morphedF0(min(iFloor+1,nAxis1))-morphedF0(iFloor));
        morphedF0onMorph(ii) = morphedF0(iFloor)+dF0;
    end;
end;
mObject3.F0 = morphedF0onMorph; 
mObject3.aperiodicityIndex = morphedApOnMorph;
mObject3.spectrogram = morphedSgramOnMorph;
mObject3.anchorTimeLocation = timeAnchorMorph(2:end-1);
%mObject3.anchorFrequency = (1-mRate)*mObject1.anchorFrequency+mRate*mObject2.anchorFrequency;
%mObject3 = fmapFromMorphto1OnTimeMorph; % This line is a dummy.
return;

%%% ------ Internal function to check for object's similarity
function mObject3 = checkForSimilarity(mObject1,mObject2)
mObject3 = [];
if mObject1.samplingFrequency ~= mObject2.samplingFrequency;mObject3 = [];return;end;
if mObject1.frameUpdateInterval ~= mObject2.frameUpdateInterval;mObject3 = [];return;end;
if length(mObject1.anchorTimeLocation) ~= length(mObject2.anchorTimeLocation);mObject3 = [];return;end;
nAnchor = length(mObject1.anchorTimeLocation);
for ii=1:nAnchor % check for similarity of anchor structure
    frequencyAnchor1 = mObject1.anchorFrequency(ii,:)';% 01/Oct./2005 by HK
    frequencyAnchor2 = mObject2.anchorFrequency(ii,:)';% 01/Oct./2005 by HK
    if (sum(frequencyAnchor1>0) ~= sum(frequencyAnchor2>0)) | ...
            (sum(frequencyAnchor1<0) ~= sum(frequencyAnchor2<0))
        display('Warning!! Object structures are inconsistent!'); % 01/Oct./2005 by HK
        return;
    end;
end;
mObject3 = createMobject;
m0bject3.samplingFrequency = mObject1.samplingFrequency;
m0bject3.frameUpdateInterval = mObject1.frameUpdateInterval;
return;

%%%--------
function mixRate = checkForMorphingConditions(mRate);
%   04/Oct./2005 added by HK

if ~isstruct(mRate)
    mixRate.F0 = mRate;
    mixRate.spectrum = mRate;
    mixRate.aperiodicity = mRate;
    mixRate.coordinate = mRate;
    return;
end;
mixRate.F0 = mRate.F0;
mixRate.spectrum = mRate.spectrum;
mixRate.aperiodicity = mRate.aperiodicity;
mixRate.coordinate = mRate.coordinate;
return;

%%%--------
function fmapFrom1to2OnTime1 = generateFrequencyMap(mObject1,mObject2);

dtFrame = mObject1.frameUpdateInterval;
endLocation1 = (length(mObject1.F0)-1)*dtFrame; % in ms
timeAnchor1 = [0;mObject1.anchorTimeLocation;endLocation1];
locationOn1 = (0:length(mObject1.F0)-1)*dtFrame;
fs = mObject1.samplingFrequency;
[nr1,nc1] = size(mObject1.spectrogram);
nAnchor = length(mObject1.anchorTimeLocation);
fmapFrom1to2 = zeros(nr1,nAnchor);
frequencyAxis = (0:nr1-1)'/(nr1-1)*fs/2;
numberOfFrequencyAnchors = zeros(nAnchor,1);
for ii=1:nAnchor
    frequencyAnchor1 = mObject1.anchorFrequency(ii,:)';
    frequencyAnchor1 = [0;frequencyAnchor1(frequencyAnchor1>0);fs/2];
    numberOfFrequencyAnchors(ii) = length(frequencyAnchor1(frequencyAnchor1>0));
    frequencyAnchor2 = mObject2.anchorFrequency(ii,:)';
    frequencyAnchor2 = [0;frequencyAnchor2(frequencyAnchor2>0);fs/2];
    fmapFrom1to2(:,ii) = interp1(frequencyAnchor1,frequencyAnchor2,frequencyAxis);
end;
for ii=1:nAnchor
    if numberOfFrequencyAnchors(ii) == 1
        if numberOfFrequencyAnchors(min(ii+1,nAnchor)) > 1
            fmapFrom1to2(:,ii) = fmapFrom1to2(:,min(ii+1,nAnchor));
        elseif numberOfFrequencyAnchors(max(ii-1,1)) > 1
            fmapFrom1to2(:,ii) = fmapFrom1to2(:,max(ii-1,1));
        end;
    end;
end;
fmapFrom1to2 = [fmapFrom1to2(:,1) fmapFrom1to2 fmapFrom1to2(:,nAnchor)];
fmapFrom1to2OnTime1 = interp1(timeAnchor1,fmapFrom1to2',locationOn1)';
return;

%%%-------
function cleanedUpObject = checkForIntegrity(inputObject);

maximumIndex = max([length(inputObject.F0), ...
    size(inputObject.spectrogram,2) ...
    size(inputObject.aperiodicityIndex,2)]);
if length(inputObject.F0) < maximumIndex
    inputObject.F0 = [inputObject.F0(:);inputObject.F0(end)*ones(maximumIndex - length(inputObject.F0),1)];
end;
if size(inputObject.spectrogram,2) < maximumIndex
    numberOfFillIn = maximumIndex-size(inputObject.spectrogram,2);
    inputObject.spectrogram = [inputObject.spectrogram inputObject.spectrogram(:,end)*ones(1,numberOfFillIn)];
end;
if size(inputObject.aperiodicityIndex,2) < maximumIndex
    numberOfFillIn = maximumIndex-size(inputObject.aperiodicityIndex,2);
    inputObject.aperiodicityIndex = [inputObject.aperiodicityIndex inputObject.aperiodicityIndex(:,end)*ones(1,numberOfFillIn)];
end;
cleanedUpObject = inputObject;

