data = load('../DATA/B18_EEG_data/EEGSleepStateData.mat');
data = data.EEGSleepStateData;

array = table2array(data(:, [1, 3]));

time = array(:,1);
state = array(:,2);

figure(1);
stairs(time, state)
xlabel('Time (mins)')
ylabel('Sleep state')
ylim([0.5, 6.5])

%%

detrend_state = state - mean(state);
z_state = detrend_state/std(detrend_state);

fs = 60;

L = length(z_state);

W = hamming(floor(L));
[ pxx, f ] = periodogram(z_state, W, [], fs);

figure(2)
plot(f, pxx)
xlabel('Frequency (CyclesHr^{-1})')
ylabel('Power Spectral Density (state^{2}(Hr)^{-1})')


