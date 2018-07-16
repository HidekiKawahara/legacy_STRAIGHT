%% Test script for optimum mixing
% 17/Feb./2016
% Hideki Kawahara

clear all
close all

nSignal = 6;
baseValue = 100;
nData = 10000;
nIteration = 5000;

stdBest = zeros(nIteration, 1);
stdMean = zeros(nIteration, 1);
stdSD = zeros(nIteration, 1);
for kk = 1:nIteration
    sdList = 0.2 + 3 * rand(nSignal, 1);
    sampleData = randn(nData, nSignal) * diag(sdList) + baseValue;
    
    H = ones(nSignal, nSignal) * baseValue^2 * 2;
    for ii = 1:nSignal
        H(ii, ii) = H(ii, ii) + 2 * sdList(ii)^2;
    end;
    v = ones(nSignal, 1) * 2 * baseValue^2;
    a = inv(H) * v;
    stdBest(kk) = std(sampleData * a);
    stdMean(kk) = std(sampleData  *ones(nSignal, 1) / nSignal);
    wSD = (1.0 ./ sdList);
    wSD = wSD / sum(wSD);
    stdSD(kk) = std(sampleData * wSD);
end;
figure;plot(stdBest, stdMean, '.');grid on;
hold all
plot([0 100],[0 100]);
axis([0 max([stdBest; stdMean]) 0 max([stdBest; stdMean])]);
axis('square');
xlabel('optimized SD');
ylabel('SD of simple mean');
figure;plot(stdBest, stdSD, '.');grid on;
hold all
plot([0 100],[0 100]);
axis([0 max([stdBest; stdSD]) 0 max([stdBest; stdSD])]);
axis('square');
xlabel('optimized SD');
ylabel('SD of 1/SD mixing');

