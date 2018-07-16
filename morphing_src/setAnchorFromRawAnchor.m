function mObject = setAnchorFromRawAnchor(mObject,rawAnchor);
%   Set anchor points using raw anchor information
%   mObject = setAnchorFromRawAnchor(mObject,rawAnchor);

%   Designed and coded by Hideki Kawahara
%   27/Feb./2005
%   Copyright(c) 2005, Hideki Kawahara

TIMINGmARGIN = 10; % threshould for merging location
[dm1,indsrt] = sort(rawAnchor(:,1));
sortedAnchor = rawAnchor(indsrt,1);
sortedFrequency = rawAnchor(indsrt,2);
indexNumber = 1:length(sortedAnchor);

%anchorCandidate = sortedAnchor(diff([-100;sortedAnchor])>TIMINGmARGIN);
anchorIndex = indexNumber(diff([-100;sortedAnchor])>TIMINGmARGIN);
anchorCandidate = sortedAnchor(anchorIndex);
mObject.anchorTimeLocation = anchorCandidate; 
nFrequency = mObject.maximumFrequencyPoints; 

nAnchor = length(anchorCandidate);
sortedAnchor(end+1) = sortedAnchor(end)+1;
anchorIndex(end+1) = anchorIndex(end)+1; % Terminator
frequencyAnchor = zeros(nAnchor,nFrequency);
for ii=1:nAnchor
    iFrequency = 0;
    anchorLocation = 0;
    for jj=1:min(nFrequency,anchorIndex(ii+1)-anchorIndex(ii)+1)
        if sortedAnchor((jj-1)+anchorIndex(ii)) < sortedAnchor(anchorIndex(ii+1))
            frequencyAnchor(ii,jj) = sortedFrequency((jj-1)+anchorIndex(ii));
            anchorLocation = anchorLocation+sortedAnchor((jj-1)+anchorIndex(ii));
            iFrequency = iFrequency+1;
        end;
    end;
    if iFrequency>1
        [dmy1,indsrt] = sort(frequencyAnchor(ii,1:iFrequency));
        frequencyAnchor(ii,1:iFrequency) = frequencyAnchor(ii,indsrt);
        mObject.anchorTimeLocation(ii) = anchorLocation/iFrequency;
    end;
end;
mObject.anchorFrequency = frequencyAnchor;
