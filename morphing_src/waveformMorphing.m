function mObject3 = waveformMorphing(mObject1,mObject2,mRate);
%   Morphing with minimum information
%   (Actually this is not real morphing.
%   It is simply blending two waveform.)
%   mObject3 = waveformMorphing(mObject1,mObject2,mRate);

%   Designed and coded by Hideki Kawahara
%   27/Feb./2005
%   Copyright(c) 2005, Hideki Kawahara

nLength = max(length(mObject1.waveform),length(mObject2.waveform));
if mObject1.samplingFrequency ~= mObject2.samplingFrequency
    mObject3 = [];
    return
end;
x = zeros(nLength,1);
x(1:length(mObject1.waveform)) = (1-mRate)*mObject1.waveform;
x(1:length(mObject2.waveform)) = mRate*mObject2.waveform + x(1:length(mObject2.waveform));

mObject3=createMobject;
mObject3.waveform = x;
mObject3.samplingFrequency = mObject1.samplingFrequency;
