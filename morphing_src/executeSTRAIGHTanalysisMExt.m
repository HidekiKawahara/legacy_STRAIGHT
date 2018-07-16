function mObject = executeSTRAIGHTanalysisMExt(mObject,optionalParameters);
%   STRAIGHT analysis for mObject
%   mObject = executeSTRAIGHTanalysisMExt(mObject,optionalParameters);
%

%   Designed and coded by Hideki Kawahara
%   26/Feb./2005
%   Copyright(c) 2005, Hideki Kawahara
%   20/March/2006 bug fix by T. Takahashi and Kawahara
%   16/Aug./2008 extended for use MulticueF0 as default

x = mObject.waveform;
fs = mObject.samplingFrequency;
if nargin>1
    %[f0raw,ap,prmF0] = exstraightsource(x,fs,optionalParameters);
    [f0raw,vuv,auxouts,prmF0]=MulticueF0v14(x,fs,optionalParameters);
    [n3sgram,analysisParamsSp]=exstraightspec(x,f0raw,fs,optionalParameters);
    [ap,analysisParams]=exstraightAPind(x,fs,f0raw,optionalParameters);
else
    %[f0raw,ap,prmF0] = exstraightsource(x,fs);
    [f0raw,vuv,auxouts,prmF0]=MulticueF0v14(x,fs);
    [n3sgram,analysisParamsSp]=exstraightspec(x,f0raw,fs);
    [ap,analysisParams]=exstraightAPind(x,fs,f0raw);
end;
if exist('vuv') % reserved for extension
    mObject.vuv = vuv;
else
    mObject.vuv = (f0raw ~= 0);
end;
temporalIndexLength=min([length(f0raw),size(n3sgram,2),size(ap,2),length(mObject.vuv)]);
mObject.F0 = f0raw(1:temporalIndexLength);
mObject.spectrogram = n3sgram(:,1:temporalIndexLength);
mObject.aperiodicityIndex = ap(:,1:temporalIndexLength);
mObject.vuv = mObject.vuv(1:temporalIndexLength);
mObject.frameUpdateInterval = prmF0.F0frameUpdateInterval;
mObject.F0extractionConditions = prmF0;
mObject.SpectrumExtractionConditions = analysisParamsSp;
mObject.AperiodicityAnalysisParams = analysisParams;
