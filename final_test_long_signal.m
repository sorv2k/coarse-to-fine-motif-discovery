% final_test_long_signal.m
% Compares brute force (gold standard) vs adaptive coarse-to-fine search
% on 200,000 point synthetic signals.
%
% Required files in same folder:
%   mpx.m
%   final_FindRescaledMotif_Adaptive.m
%   final_generate_signals_with_scale.m

clear; clc;

set(groot, 'DefaultFigureColor', 'white');
set(groot, 'DefaultAxesColor', 'white');
set(groot, 'DefaultTextColor', 'black');
set(groot, 'DefaultAxesXColor', 'black');
set(groot, 'DefaultAxesYColor', 'black');

%% ---- CONFIGURE THIS ----
signal_length = 200000;
motif_length  = 512;
true_scale    = 0.8;
rng(42);

%% ---- Generate Signals ----
fprintf('Generating synthetic signals...\n');
fprintf('Signal length: %d points\n', signal_length);
fprintf('Motif length:  %d\n', motif_length);
fprintf('True scale:    %.2f\n\n', true_scale);

[A, B] = final_generate_signals_with_scale(true_scale, signal_length, motif_length);

fprintf('Signal A: %d points\n', length(A));
fprintf('Signal B: %d points\n\n', length(B));

%% ---- BRUTE FORCE (Gold Standard) ----
fprintf('Running BRUTE FORCE (exhaustive, 101 tests at 0.01 step)...\n');
fprintf('Scale range: 0.5 to 1.5, step 0.01\n');
fprintf('Estimated time: ~%d minutes based on 50k calibration\n\n', round(1086.5*16/60));

tic;
result_bf = brute_force_search(A, B, motif_length, 0.5, 1.5, 0.01);
time_bf = toc;

fprintf('\nBrute force done.\n');
fprintf('Time:        %.1f seconds (%.2f minutes)\n', time_bf, time_bf/60);
fprintf('Found scale: %.3f\n', result_bf.bestScale);
fprintf('True scale:  %.3f\n', true_scale);
fprintf('Error:       %.2f%%\n\n', abs(result_bf.bestScale - true_scale)/true_scale*100);

%% ---- ADAPTIVE SEARCH ----
fprintf('Running ADAPTIVE COARSE-TO-FINE search...\n');
fprintf('Coarse: 11 tests (0.1 step), Fine: ~31 tests (0.01 step)\n\n');

tic;
result_ad = final_FindRescaledMotif_Adaptive(A, B, motif_length, 0.5, 1.5, 0.1, false);
time_ad = toc;

fprintf('Adaptive done.\n');
fprintf('Time:        %.1f seconds (%.2f minutes)\n', time_ad, time_ad/60);
fprintf('Found scale: %.3f\n', result_ad.bestScale);
fprintf('True scale:  %.3f\n', true_scale);
fprintf('Error:       %.2f%%\n\n', abs(result_ad.bestScale - true_scale)/true_scale*100);

%% ---- Summary ----
speedup = time_bf / time_ad;

fprintf('==========================================\n');
fprintf('           RESULTS SUMMARY\n');
fprintf('==========================================\n');
fprintf('Signal length:     %d points\n', signal_length);
fprintf('True scale factor: %.3f\n', true_scale);
fprintf('------------------------------------------\n');
fprintf('BRUTE FORCE (gold standard)\n');
fprintf('  Time:            %.1f sec (%.2f min)\n', time_bf, time_bf/60);
fprintf('  Found scale:     %.3f\n', result_bf.bestScale);
fprintf('  Error:           %.2f%%\n', abs(result_bf.bestScale-true_scale)/true_scale*100);
fprintf('------------------------------------------\n');
fprintf('ADAPTIVE SEARCH\n');
fprintf('  Time:            %.1f sec (%.2f min)\n', time_ad, time_ad/60);
fprintf('  Found scale:     %.3f\n', result_ad.bestScale);
fprintf('  Error:           %.2f%%\n', abs(result_ad.bestScale-true_scale)/true_scale*100);
fprintf('------------------------------------------\n');
fprintf('SPEEDUP:           %.2fx faster\n', speedup);
fprintf('TIME SAVED:        %.1f sec (%.2f min)\n', time_bf-time_ad, (time_bf-time_ad)/60);
fprintf('==========================================\n\n');

%% ---- Figure: Distance Curves ----
figure('Color', 'white');

plot(result_bf.scales, result_bf.distances, 'b-', 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Brute Force (%.0f sec)', time_bf));
hold on;

plot(result_ad.scales, result_ad.distances, 'r--', 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Adaptive (%.0f sec)', time_ad));

[~, bf_idx] = min(result_bf.distances);
plot(result_bf.scales(bf_idx), result_bf.distances(bf_idx), 'bs', ...
     'MarkerSize', 8, 'LineWidth', 2, ...
     'DisplayName', sprintf('BF found: %.3f', result_bf.bestScale));

[~, ad_idx] = min(result_ad.distances);
plot(result_ad.scales(ad_idx), result_ad.distances(ad_idx), 'ro', ...
     'MarkerSize', 8, 'LineWidth', 2, ...
     'DisplayName', sprintf('Adaptive found: %.3f', result_ad.bestScale));

xline(true_scale, 'g--', 'LineWidth', 1.5);
text(true_scale + 0.01, max(result_bf.distances)*0.95, ...
     sprintf('True: %.2f', true_scale), 'Color', [0 0.55 0], 'FontSize', 10);

hold off;
box off;
legend('Location', 'eastoutside');
xlabel('Scale Factor');
ylabel('Distance');
title(sprintf('Brute Force vs Adaptive  |  N=%d  |  Speedup: %.1fx', ...
              signal_length, speedup));

%% ======================================================
%% LOCAL FUNCTION: Brute Force Search with Progress Bar
%% ======================================================
function result = brute_force_search(A, B, m, MinSF, MaxSF, StepS)
    scales    = MinSF:StepS:MaxSF;
    distances = inf(numel(scales), 1);
    minlag    = round(m/2);
    split     = length(A);
    n_scales  = numel(scales);

    A = A(:);
    B = B(:);

    tic_total = tic;

    for i = 1:n_scales
        sf = scales(i);

        new_len  = max(2, round(length(B) * sf));
        B_scaled = interp1(1:length(B), B, ...
                   linspace(1, length(B), new_len), 'linear')';

        concat = [A; B_scaled];
        [mp, mpi] = mpx(concat, minlag, m);

        bestDist = inf;
        for j = 1:min(split, length(mp))
            if j + m - 1 > split, continue; end
            p = mpi(j);
            if ~isnan(p) && p > split
                bStart = p - split;
                if bStart + m - 1 > length(B_scaled), continue; end
                if std(concat(j:j+m-1)) < 1e-6, continue; end
                if mp(j) < bestDist
                    bestDist = mp(j);
                end
            end
        end

        distances(i) = bestDist;
        clear B_scaled concat mp mpi;

        % Progress bar
        elapsed   = toc(tic_total);
        avg_time  = elapsed / i;
        remaining = avg_time * (n_scales - i);
        pct       = i / n_scales * 100;
        bar_len   = 30;
        filled    = round(bar_len * i / n_scales);
        bar_str   = [repmat('=', 1, filled), repmat('-', 1, bar_len - filled)];

        % Format remaining time nicely
        if remaining >= 3600
            eta_str = sprintf('%.1f hrs', remaining/3600);
        elseif remaining >= 60
            eta_str = sprintf('%.1f min', remaining/60);
        else
            eta_str = sprintf('%.0f sec', remaining);
        end

        fprintf('[%s] %3.0f%%  Test %3d/%d  Elapsed: %5.1fs  ETA: %s\n', ...
                bar_str, pct, i, n_scales, elapsed, eta_str);
    end

    fprintf('\n');
    [bestDist, kbest] = min(distances);
    result.scales     = scales;
    result.distances  = distances;
    result.bestScale  = scales(kbest);
    result.bestDist   = bestDist;
end