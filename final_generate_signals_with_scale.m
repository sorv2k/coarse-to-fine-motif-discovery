function [A, B] = final_generate_signals_with_scale(scale_factor, length_total, motif_len)
    % Generate A with normal motif
    one_ts = zscore(cumsum(randn(1, floor(length_total/2))));
    motif_A = sin(0.1:0.01856*1.0:9.6);
    motif_A = motif_A(1:min(motif_len, length(motif_A)));  % Safety check
    one_ts = zscore([one_ts, (motif_A * 0.5) + one_ts(end)]);
    one_ts = zscore(one_ts + randn(size(one_ts))/50);
    two_ts = zscore(cumsum(randn(1, length_total - length(one_ts))));
    A = [one_ts, two_ts + (one_ts(end) - two_ts(1))];
    
    % Generate B with scaled motif
    % For scale_factor < 1.0, B needs LONGER motif (stretched)
    % For scale_factor > 1.0, B needs SHORTER motif (compressed)
    one_ts_B = zscore(cumsum(randn(1, floor(length_total/2))));
    
    % Create stretched/compressed motif for B
    step_size = 0.01856 * scale_factor;  % NOT divided, multiplied!
    motif_B = sin(0.1:step_size:9.6);
    
    % Ensure motif_B has reasonable length
    if length(motif_B) < 100
        % If too short, extend the range
        motif_B = sin(0.1:step_size:15);
    end
    
    one_ts_B = zscore([one_ts_B, (motif_B * 0.5) + one_ts_B(end)]);
    one_ts_B = zscore(one_ts_B + randn(size(one_ts_B))/50);
    two_ts_B = zscore(cumsum(randn(1, max(100, length_total - length(one_ts_B)))));
    B = [one_ts_B, two_ts_B + (one_ts_B(end) - two_ts_B(1))];
end