function mObject=createMobject
%   Create Mobjcet for morphing
%   mObject=createMobject

%   Designed and coded by Hideki Kawahara
%   25/February/2005
%   14/October/2005 Added creator information

mObject.date = datestr(now);
mObject.pwd = pwd;
mObject.waveform = [];
mObject.samplingFrequency = 44100; % default frequency
mObject.F0 = [];
mObject.vuv = [];
mObject.spectrogram = [];
mObject.aperiodicityIndex = [];
mObject.frameUpdateInterval = 1; % default frame is 1ms
mObject.anchorTimeLocation = []; 
mObject.maximumFrequencyPoints = 9; % default max frequency anchor points
mObject.anchorFrequency = [];
mObject.F0extractionConditions = [];
mObject.SpectrumExtractionConditions = [];
mObject.creatorInformation = which('createMobject');
