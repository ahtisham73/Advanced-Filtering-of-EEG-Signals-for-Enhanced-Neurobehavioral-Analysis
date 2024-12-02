% Define the path to EEGLAB and add it to MATLAB's search path
eeglabPath = 'C:\Users\ahtis\OneDrive - National University of Sciences & Technology\Semester\SEMESTER 6\DSP\LAB\project0\eeglab2024.0';
addpath(eeglabPath);

% Load the EEG dataset
eegDataset = '825_2_PD_REST.mat';
load(eegDataset);
%%
% Extract the sampling frequency (in Hz)
Fs = EEG.srate;

% Set the analysis time window in seconds
startTime = 0;
endTime = 10;

% Convert the time window to sample indices
startIndex = round((startTime - EEG.xmin) *Fs) + 1;
endIndex = round((endTime - EEG.xmin) * Fs);

% Determine the number of samples in the analysis window
sampleCount = endIndex - startIndex + 1;

% Generate the frequency vector for the FFT
freqVector = (0:(sampleCount/2)) * (Fs / sampleCount);

% Prepare a matrix to hold the PSD values for each channel
channelPSDs = zeros(length(freqVector), EEG.nbchan);

% Calculate the PSD for each channel within the specified time window
for channelIndex = 1:EEG.nbchan
    % Extract the data for the current channel and time window
    currentData = EEG.data(channelIndex, startIndex:endIndex);
    
    % Perform the Fast Fourier Transform (FFT) on the data
    fftResult = fft(currentData);
    
    % Compute the two-sided spectrum and then the single-sided spectrum
    twoSidedSpectrum = abs(fftResult / sampleCount);
    singleSidedSpectrum = twoSidedSpectrum(1:(sampleCount/2 + 1));
    singleSidedSpectrum(2:end-1) = 2 * singleSidedSpectrum(2:end-1);
    
    % Store the PSD for the current channel
    channelPSDs(:, channelIndex) = singleSidedSpectrum;
end

% Calculate the average PSD across all channels
averagePSD = mean(channelPSDs, 2);

% Find the PSD values at 50 Hz and 60 Hz
psdAt50Hz = averagePSD(round(50 / (Fs / sampleCount)) + 1);
psdAt60Hz = averagePSD(round(60 / (Fs / sampleCount)) + 1);

% Output the PSD values at 50 Hz and 60 Hz
fprintf('Average PSD at 50 Hz: %f\n', psdAt50Hz);
fprintf('Average PSD at 60 Hz: %f\n', psdAt60Hz);

% Create a plot of the average PSD
figure;
plot(freqVector, 10 * log10(averagePSD));
title('Average Power Spectral Density (PSD) Across All Channels');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
axis([0, Fs/2, -100, 50]);
grid on;
%%
% Define notch filter parameters for 60 Hz and its harmonics
notchFreqs = 60:60:(Fs/2); % Frequencies to remove
Q = 35; % Quality factor

% Design and apply notch filters for each harmonic
for freq = notchFreqs
    [b, a] = iirnotch(freq/(Fs/2), freq/(Fs/2)/Q);
    EEG.data = filtfilt(b, a, EEG.data')'; % Apply the filter to the data
end

%%
% Compute the PSD after filtering
[psd, f] = pwelch(EEG.data', hamming(EEG.pnts), [], [], Fs);

mean_psd = mean(psd, 2);
% Plot the PSD after notch filtering
figure;
plot(f, 10*log10(mean_psd));
title('Average Power Spectral Density (PSD) Across All Channels After Notch Filtering');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
axis([0, Fs/2, -100, 200]);
grid on;

% Extract the updated EEG data after notch filtering
updatedEEGData = EEG.data;

% List available channels (1 to 67)
nbchan = 1:67;

% Select 10 channels from the updated EEG data
selectedChannels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
eegSelected = updatedEEGData(selectedChannels, :);
disp('Selected Channels:');
disp(selectedChannels);

%%
%Plot the combined EEG data for selected channels
figure;
hold on;
colors = lines(10); % Get distinct colors for each plot
for i = 1:10
    plot(eegSelected(i, :), 'Color', colors(i, :));
end
hold off;
title('Combined EEG Data for Selected 10 Channels After Notch Filter');
xlabel('Time (samples)');
ylabel('Amplitude');
legend(arrayfun(@(x) ['Channel ' num2str(selectedChannels(x))], 1:10, 'UniformOutput', false));
grid on;

%%
% Generate the frequency vector for the FFT
freqVector = (0:(sampleCount/2)) * (Fs / sampleCount);

% Define filter orders
N = 1600; % Filter order

% Define the cut-off frequencies for each filter (as a fraction of the Nyquist frequency)
Fc_delta = [0.5 4] / (Fs/2); % Delta band
Fc_theta = [4 7] / (Fs/2);   % Theta band
Fc_alpha = [8 12] / (Fs/2);   % Alpha band
Fc_sigma = [12 16] / (Fs/2);  % Sigma band
Fc_beta = [13 30] / (Fs/2);   % Beta band

% Initialize filtered EEG data matrix (3D matrix)
filteredEEGData = zeros(EEG.nbchan, sampleCount, 5);
% Filter each selected EEG channel for each frequency band
for i = 1:length(selectedChannels)
    channelIndex = selectedChannels(i);
    currentData = EEG.data(channelIndex, startIndex:endIndex);
    
    % Delta band filter (combination of low pass and high pass filters)
    b_LPF_delta = fir1(N, Fc_delta(2), 'low', hanning(N+1)); % Low pass filter
    b_HPF_delta = fir1(N, Fc_delta(1), 'high', hanning(N+1)); % High pass filter
    filtered_delta = filtfilt(b_LPF_delta, 1, filtfilt(b_HPF_delta, 1, currentData));
    
    % Theta band filter (band pass filter)
    b_BPF_theta = fir1(N, Fc_theta, 'bandpass', hanning(N+1)); % Band pass filter
    filtered_theta = filtfilt(b_BPF_theta, 1, currentData);
    
    % Alpha band filter (combination of low pass, high pass, and band pass filters)
    b_LPF_alpha = fir1(N, Fc_alpha(2), 'low', hanning(N+1)); % Low pass filter
    b_HPF_alpha = fir1(N, Fc_alpha(1), 'high', hanning(N+1)); % High pass filter
    b_BPF_alpha = fir1(N, Fc_alpha, 'bandpass', hanning(N+1)); % Band pass filter
    filtered_alpha = filtfilt(b_LPF_alpha, 1, filtfilt(b_HPF_alpha, 1, filtfilt(b_BPF_alpha, 1, currentData)));
    
    % Sigma band filter (band pass filter)
    b_BPF_sigma = fir1(N, Fc_sigma, 'bandpass', hanning(N+1)); % Band pass filter
    filtered_sigma = filtfilt(b_BPF_sigma, 1, currentData);
    
    % Beta band filter (combination of low pass, high pass, and band pass filters)
    b_LPF_beta = fir1(N, Fc_beta(2), 'low', hanning(N+1)); % Low pass filter
    b_HPF_beta = fir1(N, Fc_beta(1), 'high', hanning(N+1)); % High pass filter
    b_BPF_beta = fir1(N, Fc_beta, 'bandpass', hanning(N+1)); % Band pass filter
    filtered_beta = filtfilt(b_LPF_beta, 1, filtfilt(b_HPF_beta, 1, filtfilt(b_BPF_beta, 1, currentData)));
    
    % Store the filtered data for each frequency band in the 3D matrix
    filteredEEGData(i,:,:) = [filtered_delta; filtered_theta; filtered_alpha; filtered_sigma; filtered_beta]';
end
%%

% Define the titles for each frequency band
bandTitles = {'Delta', 'Theta', 'Alpha', 'Sigma', 'Beta'};
bandColors = {'b', 'r', 'g', 'm', 'c'}; % Colors for each band

% Channels to plot results for
channelsToPlot = [1]; % Specify which channels to plot

for channelToPlot = channelsToPlot
    % Plot the filtered EEG data for the specified channels in the time domain
    figure;
    % Loop through each frequency band
    for bandIndex = 1:5
        subplot(5, 1, bandIndex);
        plot(filteredEEGData(channelToPlot, :, bandIndex), bandColors{bandIndex});
        title([bandTitles{bandIndex} ' Filtered EEG Data for Channel ' num2str(channelToPlot)]);
        xlabel('Time (samples)');
        ylabel('Amplitude');
        grid on;
    end
%%
    % Define the band names, colors, and cutoff frequencies for plotting
bands = {'Delta', 'Theta', 'Alpha', 'Sigma', 'Beta'};
colors = {'b', 'r', 'g', 'm', 'c'};
cutoffs = {[0.5, 4], [4, 7], [8, 12], [12, 16], [13, 30]};
behavioralTraits = {
    'Delta (0.5 to 4 Hz): Associated with deep sleep and restorative sleep stages.',
    'Theta (4 to 7 Hz): Related to light sleep, relaxation, creativity, and meditative states.',
    'Alpha (8 to 12 Hz): Indicative of relaxed wakefulness, reduced stress, and enhanced learning.',
    'Sigma (12 to 16 Hz): Linked to sleep spindles during sleep, involved in memory consolidation.',
    'Beta (13 to 30 Hz): Associated with active thinking, concentration, and problem-solving.'
};

% Plot the filtered EEG data for selected channels in the frequency domain
figure;
for ch = 1:length(selectedChannels)
    channelToPlot = selectedChannels(ch);
    % Loop through each band and plot
    for i = 1:length(bands)
        subplot(5, 1, i);
        % Compute FFT
        fftData = fft(filteredEEGData(channelToPlot, :, i));
        fftData = abs(fftData / sampleCount);
        fftData = fftData(1:sampleCount/2+1);
        fftData(2:end-1) = 2 * fftData(2:end-1);
        % Plot the frequency domain data
        plot(freqVector, 10*log10(fftData), colors{i});
        title([bands{i} ' Filtered EEG Data (Freq Domain) for Channel ' num2str(channelToPlot)]);
        xlabel('Frequency (Hz)');
        ylabel('Power/F (dB/Hz)');

        grid on;

        % Add markers at the cutoff frequencies
        hold on;
        plot([cutoffs{i}(1), cutoffs{i}(1)], ylim, '--k', 'LineWidth', 1); % Start cutoff frequency
        plot([cutoffs{i}(2), cutoffs{i}(2)], ylim, '--k', 'LineWidth', 1); % End cutoff frequency
        hold off;
        
        % Print behavioral traits associated with each band
        fprintf('\nBehavioral traits for %s band:\n%s\n', bands{i}, behavioralTraits{i});
    end
    pause; % Pause to view each channel's plots
end
end