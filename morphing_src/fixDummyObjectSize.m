function dummyObject = fixDummyObjectSize(dummyObject,originalObject);

frameUpdateInterval = dummyObject.frameUpdateInterval;
endMargin = size(originalObject.spectrogram,2)*frameUpdateInterval-max(originalObject.anchorTimeLocation);
if size(dummyObject.spectrogram,2)*frameUpdateInterval < max(dummyObject.anchorTimeLocation)+endMargin
    dummyFrameSize = max(dummyObject.anchorTimeLocation)+endMargin;
    dimmyFrequencySize = size(originalObject.spectrogram,1);
    dummyObject.spectrogram = ones(dimmyFrequencySize,dummyFrameSize);
    dummyObject.aperiodicityIndex = ones(dimmyFrequencySize,dummyFrameSize);
    dummyObject.F0 = ones(1,dummyFrameSize);
    dummyObject.vuv = ones(1,dummyFrameSize);
end;
