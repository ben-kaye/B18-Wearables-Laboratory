%% load data

eyesopen = load('../DATA/B18_EEG_data/EEGeyesopen.mat').eyesopen;
eyesclosed = load('../DATA/B18_EEG_data/EEGeyesclosed.mat').eyesclosed;
 
cutoff = 512;

figure(1);
tiledlayout(1,2)
nexttile
plot(eyesopen(1:cutoff))
xlabel('sample #')
ylabel('amplitude (\muV)')
title('Eyes open')
xlim([0 512])

nexttile
plot(eyesclosed(1:512))
xlabel('sample #')
ylabel('amplitude (\muV)')
title('Eyes closed')
xlim([0 512])

%% Processing
eyesopen_detrend = standardise(eyesopen);
eyesclosed_detrend = standardise(eyesclosed);

fs = 256; % Sampling frequency [Hz]

eyesopen_detrend = standardise(eyesopen);
eyesclosed_detrend = standardise(eyesclosed);

% LPF at 40 Hz to get rid of interfence sources + noise
fc = 40; % [Hz]
Norder = 4;
[ B, A ] = butter(Norder, fc/fs);
filt_eyesopen = filter(B, A, eyesopen_detrend);
filt_eyesclosed = filter(B, A, eyesclosed_detrend);

filt_eyesopen = filt_eyesopen(1:cutoff);
filt_eyesclosed = filt_eyesclosed(1:cutoff);

figure(3);
tiledlayout(1,2)
nexttile
plot(filt_eyesopen(1:cutoff))
xlabel('sample #')
ylabel('processed amplitude')
title('Eyes open')
xlim([0 512])
nexttile
plot(filt_eyesclosed(1:512))
xlabel('sample #')
ylabel('processed amplitude')
title('Eyes closed')
xlim([0 512])


%% Spectral analysis

signals = [ filt_eyesopen; filt_eyesclosed ];
[ p_xx, f_xx ] = periodogram(signals', hamming(cutoff), [], fs, 'psd');

figure(2)
tiledlayout(1,2)
nexttile
plot(f_xx, p_xx(:,1))
xlim([0, 60])
xlabel('frequency (Hz)')
ylabel('power spectral density \muV^{2}Hz^{-1})')
title('Eyes open')
% ylim([0, 5e-2])
nexttile
plot(f_xx, p_xx(:,2))
xlabel('frequency (Hz)')
ylabel('power spectral density \muV^{2}Hz^{-1})')
title('Eyes closed')
xlim([0, 60])
% ylim([0, 1e-2])

function standard_data = standardise(data)
    N_samples = length(data);
    data = data - mean(data);
    sd = sqrt(sum(data.^2)/(N_samples - 1));
    standard_data = data/sd;
end
