function syntheszed_signal = SynthesizeLegacy_STRAIGHT_default(x, fs)
% Conditions are based on the web document
% 

f0raw = MulticueF0v14(x,fs);
ap = exstraightAPind(x,fs,f0raw);
n3sgram=exstraightspec(x,f0raw,fs);
syntheszed_signal = exstraightsynth(f0raw,n3sgram,ap,fs);
end