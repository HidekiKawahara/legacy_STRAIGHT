function displayMobject(mObject,fieldname,note)
%   M-object information display
%   displayMobject(mObject,fieldname,note);

%   Designed and coded by Hideki Kawahara
%   27/Feb./2005
%   Copyright(c) 2005, Hideki Kawahara
%   05/Oct./2005 minor bug fix

fs = mObject.samplingFrequency;
tFrame = mObject.frameUpdateInterval;
switch fieldname
    case 'spectrogram'
        figure
        [nrow,ncolumn]=size(mObject.spectrogram);
        timeSpan = [0 (ncolumn-1)*tFrame];
        dBsgram = 20*log10(mObject.spectrogram);
        maxSgramdB = max(max(dBsgram));
        imagesc(timeSpan, [0 fs/2],max(dBsgram,maxSgramdB-70));
        axis('xy');
        set(gca,'fontsize',14);
        xlabel('time (ms)');
        ylabel('frequency (Hz)');
        title([note ' time span 0 ' num2str(timeSpan(2),10) ' (ms) ' datestr(now)]);
    case 'waveform'
        figure
        x = mObject.waveform;
        timeSpan = (0:length(x)-1)/fs*1000;
        plot(timeSpan,x);grid on;
        axis([timeSpan(1) timeSpan(end) 1.1*[min(x) max(x)]]);
        set(gca,'fontsize',14);
        xlabel('time (ms)');
        title([note ' time span 0 ' num2str(round(timeSpan(end)),8) ' (ms) ' datestr(now)]);
    case {'anchorFrequency', 'anchorTimeLocation'}
        figure
        [nrow,ncolumn]=size(mObject.spectrogram);
        timeSpan = [0 (ncolumn-1)*tFrame];
        dBsgram = 20*log10(mObject.spectrogram);
        maxSgramdB = max(max(dBsgram));
        imagesc(timeSpan, [0 fs/2],max(dBsgram,maxSgramdB-70));
        axis('xy');
        set(gca,'fontsize',14);
        xlabel('time (ms)');
        ylabel('frequency (Hz)');
        title([note ' time span 0 ' num2str(timeSpan(2),10) ' (ms) ' datestr(now)]);
        if length(mObject.anchorTimeLocation)>0
            hold on;
            for ii=1:length(mObject.anchorTimeLocation)
                hh = plot(mObject.anchorTimeLocation(ii)*[1 1],[0 fs/2],'w:');
                set(hh,'linewidth',2);
                if sum(mObject.anchorFrequency(ii,:)>0)>0
                    nFrequency = sum(mObject.anchorFrequency(ii,:)>0);
                    anchorFrequencyVector = mObject.anchorFrequency(ii,mObject.anchorFrequency(ii,:)>0); % 05/Oct./2005 HK
                    for jj=1:nFrequency
                        hh=plot(mObject.anchorTimeLocation(ii),anchorFrequencyVector(jj),'ok');
                        set(hh,'markersize',9);
                        set(hh,'linewidth',2);
                        hh=plot(mObject.anchorTimeLocation(ii),anchorFrequencyVector(jj),'.w');
                        set(hh,'markersize',7);
                        set(hh,'linewidth',4);
                    end;
                end;
            end;
            hold off;
        end;
end;

