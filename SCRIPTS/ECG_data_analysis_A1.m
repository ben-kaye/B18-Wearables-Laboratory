%% Initialise data

ECG_data = load('../DATA/B18_ECG_data/PhysionetData.mat');
ECG_metadata = ECG_data.RecordingInfo;
ECG_signals = ECG_data.Signals;

%% Display table data

figure(1)
histogram(ECG_metadata.Label)
ylabel('Count')

%% Ex1.1

% using A05498 as normal subject :: idx 3443
normal_subject_index = 3443;
ECG_cell = ECG_signals(normal_subject_index); 
ECG_N = 1e-2*ECG_cell{1}; % normalise

fs = 300; % sample rate {Hz}

% locations P, Q, R, S, T x2 
xy_idx = [ 3008, 3022, 3044, 3052, 3127, 3414, 3331, 3339, 3312, 3297 ];
xy_pqrst = [ 10.03, 0.16; 10.07, -0.39; 10.15, 3.31; 10.17, -2.03; ...
    10.42, 2.30; 10.99, 0.22; 11.04, -0.23; 11.10, 3.34; 11.13, -2.55; 11.38, 2.34 ];



sample_length = size(ECG_N,2);

time = (1:size(ECG_N,2))/fs; % generate time

plot(time, ECG_N, '-o', 'MarkerIndices', xy_idx, 'MarkerEdgeColor', 'r')
xlim([9.7, 11.7]);
xlabel('Time (s)')
ylabel('Amplitude (mV)')
text(xy_pqrst(:, 1) + 0.02, xy_pqrst(:,2) - 0.1, {'P', 'Q', 'R', 'S', 'T', 'P', 'Q', 'R', 'S', 'T'}, 'Color','#E4951B')

%% Ex1.2


% extract data
subject_indices = [ 3443, 3016, 3017, 3142 ];
cat = ECG_metadata{subject_indices, 2};
subject_id = ECG_metadata{subject_indices, 1};
titles = {'Normal rhythmm', 'Abnormal-fibrillation', 'Other rhythm', 'Noisy'};
labels = ECG_metadata{subject_indices, 1};
ECG_traces = 1e-2*[ECG_signals{subject_indices}]; % include normalisation factor
ECG_traces = reshape(ECG_traces, [ sample_length, 4 ])';



figure(2)
Tplot = tiledlayout(4,1);
for k_plot = 1:4
    nexttile
    plot(time, ECG_traces(k_plot, :))
    xlabel('Time (s)')
    ylabel('Amplitude (mV)')
    title(sprintf('Type %s; Recording %s', titles{k_plot}, subject_id{k_plot})); 
    
end