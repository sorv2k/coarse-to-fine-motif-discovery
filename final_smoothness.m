%% Curve Smoothness Analysis
clear; clc;

fprintf('=== DISTANCE CURVE SMOOTHNESS ANALYSIS ===\n\n');

test_scales = [0.6, 0.75, 0.8, 1.0, 1.2, 1.4];

figure('Color', 'white');

for i = 1:length(test_scales)
    true_scale = test_scales(i);
    
    rng(1000 + i);
    [A, B] = final_generate_signals_with_scale(true_scale, 20000, 512);
    
    % Get full distance curve
    fprintf('Computing curve for scale %.2f...\n', true_scale);
    result = final_FindRescaledMotif_Adaptive(A, B, 512, 0.5, 1.5, 0.01, false);
    
    subplot(2, 3, i);
    plot(result.scales, result.distances, 'b-', 'LineWidth', 1.5);
    box off;
    xlabel('Scale Factor');
    ylabel('Distance');
    title(sprintf('True Scale = %.2f', true_scale));
    
    % Mark true scale location
    [~, true_idx] = min(abs(result.scales - true_scale));
    
    % Calculate smoothness
    dists_norm = (result.distances - min(result.distances)) / ...
                 (max(result.distances) - min(result.distances));
    smoothness = std(diff(dists_norm));
    
    end

sgtitle('Distance Curves for Different Embedded Scales', 'FontWeight', 'bold');