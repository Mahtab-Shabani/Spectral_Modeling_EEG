% EEG_Spectral_Analysis_2014a_compatible.m
% MATLAB R2014a compatible script to load EEG.mat and compute smoothed spectrum + peak detection

clear; close all; clc;

%% ========== USER SETTINGS ==========
matfile = 'EEG.mat';     % put EEG.mat in current folder or give full path
Fs = 250;                % SAMPLING FREQUENCY in Hz (CHANGE to real value if known)
smooth_win = 7;          % smoothing window (odd integer, samples in frequency domain)
peak_threshold_std = 2;  % threshold = mean(spec_smooth) + peak_threshold_std * std(spec_smooth)
use_signal_name = '';    % if non-empty, will try to use a variable by that name (e.g., 'alpha')
channel_index = 1;       % if using raw multi-channel data, which channel to use (1-based)
plot_xlim = [0 60];      % frequency range to show in plots (Hz) - adjust as needed

%% ========== LOAD DATA ==========
S = load(matfile);
varnames = fieldnames(S);
fprintf('Variables in %s:\n', matfile);
for i=1:length(varnames)
    v = varnames{i};
    sz = size(S.(v));
    fprintf('  %s    size=%s\n', v, mat2str(sz));
end

% Decide signal to analyze:
% Priority: user-specified name -> 'alpha' vector (if exists) -> raw matrix x12120x2Eedf column -> first 5000-length band delta
signal = [];
if ~isempty(use_signal_name) && isfield(S,use_signal_name)
    signal = S.(use_signal_name);
    fprintf('Using user-specified variable: %s\n', use_signal_name);
elseif isfield(S,'alpha')
    signal = S.alpha;
    fprintf('Using variable: alpha\n');
elseif isfield(S,'x12120x2Eedf')
    raw = S.x12120x2Eedf;
    if size(raw,2) >= channel_index
        signal = raw(:,channel_index);   % choose channel
        fprintf('Using raw multichannel: x12120x2Eedf, channel %d\n', channel_index);
    else
        error('Requested channel_index exceeds number of channels in raw data.');
    end
elseif isfield(S,'dalt')
    signal = S.dalt;
    fprintf('Using variable: dalt (delta band)\n');
else
    error('No suitable signal variable found in the .mat file. Set use_signal_name to one of the variable names shown above.');
end

% Ensure column vector
signal = signal(:);
N = length(signal);
fprintf('Selected signal length: %d samples\n', N);

%% ========== OPTIONAL: trim / resample if needed ==========
% If your Fs assumption does not match data length/expected time, change Fs above.
% You can also downsample/resample here if desired.

%% ========== TIME-DOMAIN PLOT ==========
t = (0:N-1)/Fs;
figure;
plot(t, signal);
xlabel('Time (s)'); ylabel('Amplitude'); title('Time-domain signal (full)');
xlim([t(1) t(end)]);

%% ========== SPECTRUM: FFT -> single-sided amplitude ==========
% Zero-mean signal improves spectral estimate
sig0 = signal - mean(signal);

% Choose FFT length (next power of two for nicer freq resolution — optional)
L = 2^nextpow2(N);
Y = fft(sig0, L);
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:floor(L/2))/L;

% Plot raw spectrum
figure;
plot(f, P1);
xlabel('Frequency (Hz)'); ylabel('|P1(f)|'); title('Raw single-sided amplitude spectrum');
xlim(plot_xlim);

%% ========== SMOOTH SPECTRUM (moving average via conv) ==========
w = ones(1, smooth_win)/smooth_win;
spec_smooth = conv(P1, w, 'same');  % R2014a-compatible

figure;
plot(f, P1, 'Color',[0.7 0.7 0.7]); hold on;
plot(f, spec_smooth, 'b', 'LineWidth', 1.2);
xlabel('Frequency (Hz)'); ylabel('Amplitude');
legend('Raw','Smoothed'); title('Smoothed Spectrum');
xlim(plot_xlim);

%% ========== PEAK DETECTION (2014a-safe) ==========
% findpeaks on the smoothed spectrum (single input)
[pks, locs] = findpeaks(spec_smooth);   % pks: amplitudes, locs: indices into spec_smooth

% convert indices -> frequency
peak_freqs = f(locs);

% compute manual threshold and filter peaks
threshold = mean(spec_smooth) + peak_threshold_std * std(spec_smooth);
mask = pks >= threshold;
peak_vals = pks(mask);
peak_freqs = peak_freqs(mask);

% if no peaks pass threshold, relax threshold automatically (optional)
if isempty(peak_vals)
    fprintf('No peaks passed threshold (%.3e). Relaxing threshold to mean + 1*std.\n', threshold);
    threshold = mean(spec_smooth) + 1*std(spec_smooth);
    mask = pks >= threshold;
    peak_vals = pks(mask);
    peak_freqs = peak_freqs(mask);
end

% Plot detected peaks
hold on;
plot(peak_freqs, peak_vals, 'ro', 'MarkerFaceColor','r');
for i=1:length(peak_freqs)
    text(peak_freqs(i)+0.5, peak_vals(i), sprintf('%.2f Hz', peak_freqs(i)));
end

% Print peaks
fprintf('Detected peaks (frequency [Hz] : amplitude):\n');
for i=1:length(peak_freqs)
    fprintf('  %.3f Hz : %.6e\n', peak_freqs(i), peak_vals(i));
end

%% ========== OPTIONAL: SAVE RESULTS ==========
% save('EEG_spectral_results.mat', 'Fs', 'N', 'f', 'P1', 'spec_smooth', 'peak_freqs', 'peak_vals', 'threshold');
% fprintf('Results saved to EEG_spectral_results.mat\n');
