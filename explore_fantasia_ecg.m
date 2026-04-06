% explore_fantasia_ecg.m
% Extracts ECG signals for one young and one elderly subject from the
% Fantasia CSV dataset, then runs FindRescaledMotif to detect temporal scaling.
%
% Dataset: fantasia_ecg_respiration_signals.csv (3.95 GB)
% Columns: ECG, RESP, Participant, Sample, Sampling_Rate, Database
% Participants: Fantasia_f1y01..f1y10 (young), Fantasia_f1o01..f1o10 (elderly)
% Sampling rate: 250 Hz

clear; clc;

set(groot, 'DefaultFigureColor', 'white');
set(groot, 'DefaultAxesColor', 'white');
set(groot, 'DefaultTextColor', 'black');
set(groot, 'DefaultAxesXColor', 'black');
set(groot, 'DefaultAxesYColor', 'black');

%% ---- CONFIGURE THIS ----
csv_path      = '/Users/souravg/Documents/MATLAB/Rastamat/Cricket_Explore/fantasia_ecg_respiration_signals.csv';
young_subject = 'Fantasia_f1y01';
old_subject   = 'Fantasia_f1o01';
N             = 20000;

%% ---- STEP 1: Extract ECG for both subjects ----
fprintf('Reading CSV for subject: %s\n', young_subject);
fprintf('(This may take 1-3 minutes for a 4GB file...)\n');
ecg_young = extract_ecg(csv_path, young_subject, N);

fprintf('Reading CSV for subject: %s\n', old_subject);
ecg_old = extract_ecg(csv_path, old_subject, N);

fprintf('Done reading.\n');
fprintf('Young subject ECG length: %d points\n', length(ecg_young));
fprintf('Old   subject ECG length: %d points\n', length(ecg_old));

%% ---- STEP 2: Estimate heart rates ----
[~, locs_young] = findpeaks(ecg_young, 'MinPeakDistance', 150, 'MinPeakProminence', 0.5);
[~, locs_old]   = findpeaks(ecg_old,   'MinPeakDistance', 150, 'MinPeakProminence', 0.5);

duration_sec = N / 250;
hr_young     = (length(locs_young) / duration_sec) * 60;
hr_old       = (length(locs_old)   / duration_sec) * 60;
expected_sf  = hr_old / hr_young;

fprintf('\n=== Heart Rate Estimates ===\n');
fprintf('Young subject (%s): %.1f bpm\n', young_subject, hr_young);
fprintf('Old   subject (%s): %.1f bpm\n', old_subject,   hr_old);
fprintf('Expected scaling factor (old/young): %.3f\n', expected_sf);

if length(locs_young) > 1
    avg_beat_len = round(mean(diff(locs_young)));
else
    avg_beat_len = round(250 * 60 / hr_young);
end
fprintf('Estimated beat length (young): %d samples\n', avg_beat_len);

%% ---- FIGURE 1: Raw ECG signals (zoomed to 2000 points) ----
figure('Color', 'white');
subplot(2,1,1);
plot(ecg_young(1:2000), 'b');
box off;
title(sprintf('%s - ECG (first 8 seconds)', young_subject));
xlabel('Sample (250 Hz)');
ylabel('ECG (mV)');

subplot(2,1,2);
plot(ecg_old(1:2000), 'r');
box off;
title(sprintf('%s - ECG (first 8 seconds)', old_subject));
xlabel('Sample (250 Hz)');
ylabel('ECG (mV)');

%% ---- FIGURE 2: Zoomed view showing individual beats ----
% 800 samples = 3.2 seconds, shows ~4 young beats and ~3 elderly beats
figure('Color', 'white');
zoom_range = 1:800;

subplot(2,1,1);
plot(zscore(ecg_young(zoom_range)), 'b');
box off;
title(sprintf('Young subject %s (%.0f bpm) - 3.2 seconds', young_subject, hr_young));
xlabel('Sample');
ylabel('Z-normalized ECG');

subplot(2,1,2);
plot(zscore(ecg_old(zoom_range)), 'r');
box off;
title(sprintf('Elderly subject %s (%.0f bpm) - 3.2 seconds', old_subject, hr_old));
xlabel('Sample');
ylabel('Z-normalized ECG');

%% ---- STEP 3: Run FindRescaledMotif ----
motif_length = avg_beat_len * 3;

fprintf('\nRunning FindRescaledMotif...\n');
fprintf('Signal A: %s (%d points)\n', young_subject, length(ecg_young));
fprintf('Signal B: %s (%d points)\n', old_subject,   length(ecg_old));
fprintf('Motif length: %d, Scale range: 0.5 to 1.5\n', motif_length);

result = final_FindRescaledMotif_Adaptive(ecg_young, ecg_old, motif_length, 0.5, 1.5, 0.01, false);

fprintf('\n=== FindRescaledMotif Result ===\n');
fprintf('Best scaling factor found: %.3f\n', result.bestScale);
fprintf('Expected scaling factor:   %.3f\n', expected_sf);
fprintf('Error: %.2f%%\n', abs(result.bestScale - expected_sf) / expected_sf * 100);

%% ---- FIGURE 4: Distance curve ----
figure('Color', 'white');

plot(result.scales, result.distances, 'b-', 'LineWidth', 1.5);
box off;
xlabel('Scale Factor');
ylabel('Distance');
title(sprintf('Distance Curve: %s vs %s', young_subject, old_subject));

%% ---- HELPER FUNCTION ----
function ecg = extract_ecg(csv_path, participant_name, N)
    fid = fopen(csv_path, 'r');
    if fid == -1
        error('Cannot open file: %s', csv_path);
    end

    header = fgetl(fid);
    cols     = strsplit(header, ',');
    ecg_col  = find(strcmpi(strtrim(cols), 'ECG'));
    part_col = find(strcmpi(strtrim(cols), 'Participant'));

    if isempty(ecg_col) || isempty(part_col)
        fclose(fid);
        error('Could not find ECG or Participant columns. Header: %s', header);
    end

    ecg   = zeros(N, 1);
    count = 0;

    while ~feof(fid) && count < N
        line = fgetl(fid);
        if ~ischar(line), break; end

        parts = strsplit(line, ',');
        if length(parts) < max(ecg_col, part_col)
            continue;
        end

        if strcmp(strtrim(parts{part_col}), participant_name)
            val = str2double(parts{ecg_col});
            if ~isnan(val)
                count = count + 1;
                ecg(count) = val;
            end
        end
    end

    fclose(fid);
    ecg = ecg(1:count);

    if count < N
        fprintf('Warning: only found %d rows for %s (requested %d)\n', ...
                count, participant_name, N);
    end
end


save('fantasia_ecg_results.mat', 'ecg_young', 'ecg_old', 'result', ...
     'hr_young', 'hr_old', 'expected_sf', 'young_subject', 'old_subject');
fprintf('Results saved to fantasia_ecg_results.mat\n');