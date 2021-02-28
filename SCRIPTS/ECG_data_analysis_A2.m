%% Initialise data

ECG_data = load('../DATA/B18_ECG_data/PhysionetData.mat');
ECG_metadata = ECG_data.RecordingInfo;
ECG_signals = ECG_data.Signals;

sample_length = 9e3;

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
org_signal; % the original unfiltered ecg signal
fs = 300; % sampling frequency [Hz]
f_low = 100; % low - pass frequency [Hz]
f_high = 0.5; % high - pass frequency [Hz]
norder = 4; % filter order

% filter_ecg .m filters the orginal signal and outputs
%the band - pass filtered ECG signal
filtered_signal = filter_ecg ( ECG_normalised, fs, f_high, f_low, norder);


