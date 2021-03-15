%% Initialise data

ECG_data = load('../DATA/B18_ECG_data/PhysionetData.mat');
ECG_metadata = ECG_data.RecordingInfo;
ECG_signals = ECG_data.Signals;

sample_length = 9e3;
fs = 300;
time = (1:sample_length)/fs; % generate time


subject_indices = [ 3443, 3016, 3017, 3142 ];
cat = ECG_metadata{subject_indices, 2};
subject_id = ECG_metadata{subject_indices, 1};

titles = {'Normal rhythmm', 'Abnormal-fibrillation', 'Other rhythm', 'Noisy'};
labels = ECG_metadata{subject_indices, 1};
ECG_traces = 1e-2*[ECG_signals{subject_indices}]; % include normalisation factor
ECG_traces = reshape(ECG_traces, [ sample_length, 4 ])';


%% ExB
means = mean(ECG_traces, 2);
detrend = ECG_traces - means;
st_devs = sqrt(sum(detrend.^2, 2)/(sample_length - 1));
ECG_normalised = detrend./st_devs; % apply sample deviation detrend

% filter parameters
fs = 300; % sampling frequency [Hz]
f_low = 40; % low - pass frequency [Hz]
f_high = 10; % high - pass frequency [Hz]
norder = 4; % filter order

% filter_ecg .m filters the orginal signal and outputs
%the band - pass filtered ECG signal
filtered_signal = filter_ecg (ECG_normalised', fs, f_high, f_low, norder)';

%% ex3

[ p_xx , f_psd ] = pwelch (filtered_signal',[] ,[] ,[] , fs ,'psd');
p_xx = p_xx';


figure(2)

PSDplot = tiledlayout(1,2);

psd_titles = {'Welchs method Periodogram (normal)', 'Welchs method Periodogram (AF)'};

for k_plot = 1:2
    nexttile
    plot(f_psd, p_xx(k_plot, :))
    xlim ([0 , f_low ])
    xlabel ('Frequency (Hz)')
    ylabel ('Power Spectral Density (mV ^{2} Hz ^{ -1}) ')
    title(sprintf('%s; Recording %s', titles{k_plot}, subject_id{k_plot})); 
    
end


%% identify peaks (ex 4)

prominence = 0.8;
min_dist = 200;% corresponds to indices in 0.7s and max identifiable heart rate of 90bpm this is fine assuming resting

% [ norm_peaks, NPI ] = findpeaks(filtered_signal(1,:), 'MinPeakProminence', prominence, 'MaxPeakWidth', 2); % norm peak index (NPI)
[ norm_peaks, NPI ] = findpeaks(filtered_signal(1,:), 'MinPeakDistance', min_dist, 'MinPeakProminence', prominence); % norm peak index (NPI)'

[ af_peaks, AFPI ] = findpeaks(filtered_signal(2,:), 'MinPeakDistance', min_dist, 'MinPeakProminence', prominence);
indices = {NPI, AFPI};

RRi = cell(1,2); % init cells
RRi_detrend = RRi;
t_RRi = RRi;
pxx_plombs = RRi;
fxx_plombs = RRi;

% have to do cell based as num of heart beats could vary between subjects.

f_interest = 1e-3:1e-3:6e-1;
for ki = 1:2 
%     temp_idcs = indices{ki}
    RRi{ki} = 1000*diff(indices{ki})/fs; % get gradient and convert to ms
    
    t_RRi{ki} = time(indices{ki}(2:end));
    
%     RRi_detrend{ki} = RRi{ki} - mean(RRi{ki});
    
    [ pxx_plomb, f_plomb ] = plomb(RRi{ki}/1000, t_RRi{ki}, f_interest, 'psd');
    
    pxx_plombs{ki} = pxx_plomb;
    fxx_plombs{ki} = f_plomb;
end

%% plot

BinWidth = 20;
for k_plot = 1:2
    figure(4 + k_plot);
    TL = tiledlayout('flow');
    title(TL, sprintf('Subject %s; Type %s', subject_id{k_plot}, titles{k_plot}))
    nexttile
    plot(time, filtered_signal(k_plot, :), '-o', 'MarkerIndices', indices{k_plot}, 'MarkerEdgeColor','red')
    xlabel('Time (s)')
    ylabel('Normalised amplitude')
    ylim([-2.5, 3])
    title('Labelled ECG')
    
    nexttile
    plot(t_RRi{k_plot}, RRi{k_plot}, '--ok')
    xlabel('Time (s)')
    ylabel('R-R interval (ms)')
    ylim([500 1450])
    title('R-R intervals')
    
    nexttile
    histogram(RRi{k_plot}, 'BinWidth', BinWidth)
    xlabel('R-R interval (ms)')
    ylabel('Count')
    xlim([500 1450])
    ylim([0 9]);
    title('Histogram of R-R intervals')
    
    nexttile
    plot(fxx_plombs{k_plot}, pxx_plombs{k_plot})
    xlabel('Frequency (Hz)')
    ylabel('Power Spectral Density (s^{2}Hz^{-1})')
    title('Lomb-Scargle PSD Spectrum')
end





