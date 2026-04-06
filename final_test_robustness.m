%% Validate Across Multiple Scenarios
clear; clc;

set(groot, 'DefaultFigureColor', 'white');
set(groot, 'DefaultAxesColor', 'white');
set(groot, 'DefaultTextColor', 'black');
set(groot, 'DefaultAxesXColor', 'black');
set(groot, 'DefaultAxesYColor', 'black');

fprintf('=== ROBUSTNESS TESTING ===\n\n');

test_scales = [0.6, 0.75, 0.8, 1.0, 1.2, 1.4];
results = struct();

figure('Color', 'white');

for i = 1:length(test_scales)
    true_scale = test_scales(i);
    fprintf('--- Test %d: True Scale = %.2fx ---\n', i, true_scale);

    rng(1000 + i);
    [A, B] = final_generate_signals_with_scale(true_scale, 20000, 512);

    % Phase 1: Coarse
    tic;
    result_coarse = final_FindRescaledMotif_Adaptive(A, B, 512, 0.5, 1.5, 0.1, false);
    time_coarse = toc;

    % Phase 2: Fine
    window = 0.15;
    MinSF = max(0.5, result_coarse.bestScale - window);
    MaxSF = min(1.5, result_coarse.bestScale + window);

    tic;
    result_fine = final_FindRescaledMotif_Adaptive(A, B, 512, MinSF, MaxSF, 0.01, false);
    time_fine = toc;

    total_time = time_coarse + time_fine;
    error_pct = abs(result_fine.bestScale - true_scale) / true_scale * 100;

    fprintf('Found: %.3fx, Error: %.2f%%, Time: %.2fs\n\n', ...
            result_fine.bestScale, error_pct, total_time);

    results(i).true_scale    = true_scale;
    results(i).found_scale   = result_fine.bestScale;
    results(i).error_pct     = error_pct;
    results(i).time          = total_time;
    results(i).scales_coarse = result_coarse.scales;
    results(i).dist_coarse   = result_coarse.distances;
    results(i).scales_fine   = result_fine.scales;
    results(i).dist_fine     = result_fine.distances;

    % Plot distance curve for this test case
    subplot(2, 3, i);
    plot(result_coarse.scales, result_coarse.distances, 'b-', 'LineWidth', 1.5);
    box off;
    xlabel('Scale Factor');
    ylabel('Distance');

    if error_pct < 5
        title(sprintf('True = %.2f  |  Found = %.3f  |  PASS', ...
              true_scale, result_fine.bestScale));
    else
        title(sprintf('True = %.2f  |  Found = %.3f  |  FAIL', ...
              true_scale, result_fine.bestScale));
    end
end

sgtitle('Adaptive Search: Distance Curves for 6 Test Cases', 'FontWeight', 'bold');

%% Summary
fprintf('=== SUMMARY OF ALL TESTS ===\n');
fprintf('%-10s %-12s %-12s %-10s\n', 'True', 'Found', 'Error', 'Time(s)');
fprintf('%-10s %-12s %-12s %-10s\n', '----', '-----', '-----', '-------');
for i = 1:length(results)
    fprintf('%-10.2f %-12.3f %-12.2f%% %-10.2f\n', ...
            results(i).true_scale, results(i).found_scale, ...
            results(i).error_pct, results(i).time);
end

avg_time    = mean([results.time]);
avg_error   = mean([results.error_pct]);
success     = sum([results.error_pct] < 5);
speedup     = 170.87 / avg_time;

fprintf('\nAverage time:    %.2f seconds\n', avg_time);
fprintf('Average error:   %.2f%%\n', avg_error);
fprintf('Success rate:    %d/%d (%.0f%%)\n', success, length(results), success/length(results)*100);
fprintf('Speedup vs BF:   %.2fx\n', speedup);