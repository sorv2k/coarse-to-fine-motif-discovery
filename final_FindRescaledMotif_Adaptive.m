function result = FindRescaledMotif_Adaptive(A, B, m, MinSF, MaxSF, StepS, PlotFlag)
    % Adaptive version with configurable search parameters
    
    if nargin < 7
        PlotFlag = false;
    end
    
    if ~isvector(A) || ~isvector(B)
        error('A and B must be vectors');
    end
    if m < 4
        error('m must be >= 4');
    end

    A = A(:);
    B = B(:);

    minlag = round(m/2);

    scales = MinSF:StepS:MaxSF;
    distances = inf(numel(scales),1);
    idxA_all  = nan(numel(scales),1);
    idxB_all  = nan(numel(scales),1);

    split = length(A);

    if PlotFlag
        fprintf('Time series A is of length %d\n', length(A));
        fprintf('Time series B is of length %d\n', length(B));
        fprintf('Motif length is %d\n', m);
        fprintf('Testing scales from %.2f to %.2f (step %.3f)\n\n', ...
                MinSF, MaxSF, StepS);
    end

    % ----- brute-force scale loop -----
    for i = 1:numel(scales)
        sf = scales(i);

        if PlotFlag
            fprintf('Testing with B scaled to %.2f\n', sf);
        end

        % Scale B
        B_scaled = rescale_linear(B, sf);

        % MPX on concatenated signal
        concat = [A; B_scaled];
        [mp, mpi] = mpx(concat, minlag, m);

        bestDist = inf;
        bestA = NaN;
        bestB = NaN;

        maxStartA = min(split, length(mp));
        for j = 1:maxStartA
            if j + m - 1 > split
                continue;
            end

            p = mpi(j);
            if ~isnan(p) && p > split
                bStart = p - split;
                if bStart + m - 1 > length(B_scaled)
                    continue;
                end

                subA = concat(j:j+m-1);
                if std(subA) < 1e-6
                    continue;
                end

                if mp(j) < bestDist
                    bestDist = mp(j);
                    bestA = j;
                    bestB = bStart;
                end
            end
        end

        distances(i) = bestDist;
        idxA_all(i)  = bestA;
        idxB_all(i)  = bestB;

        clear B_scaled concat mp mpi;
    end

    % ----- best scale -----
    [bestDistFinal, kbest] = min(distances);
    bestScale = scales(kbest);
    bestIdxA  = idxA_all(kbest);
    bestIdxB_scaled = idxB_all(kbest);

    bestIdxB_unscaled = round(bestIdxB_scaled / bestScale);
    bestIdxB_unscaled = max(1, min(bestIdxB_unscaled, length(B)-m+1));

    % ----- prints -----
    if PlotFlag
        fprintf('\nThe motif is at location %d in A, and location %d in B\n', ...
                bestIdxA, bestIdxB_unscaled);
        fprintf('The best match required the data in B to be rescaled to %.0f%%\n', ...
                100*bestScale);
        fprintf('The Euclidean distance after rescaling was %.4f\n', bestDistFinal);
    end

    % ----- plots -----
    if PlotFlag
        % Plot 1: distance curve
        figure('Color','white');
        plot(scales, distances, 'b-', 'LineWidth', 2);
        xlabel('Scale Factor');
        ylabel('Distance');
        title('Distance vs Scale');

        % Plot 2: 2x1 motif visualization
        B_best = rescale_linear(B, bestScale);

        motifA = zscore(extract_window(A, bestIdxA, m));
        motifB_unscaled = zscore(extract_window(B, bestIdxB_unscaled, m));
        motifB_scaled   = zscore(extract_window(B_best, bestIdxB_scaled, m));

        dist_before = sqrt(sum((motifA - motifB_unscaled).^2));
        dist_after  = sqrt(sum((motifA - motifB_scaled).^2));

        figure('Color','white');
        sgtitle(sprintf('Motif Visualization (Found: %.3fx)', ...
                 bestScale), 'FontWeight','bold');

        subplot(2,1,1);
        plot(motifA,'b','LineWidth',2); hold on;
        plot(motifB_unscaled,'r','LineWidth',2);
        title(sprintf('BEFORE Scaling (Distance: %.2f)', dist_before), ...
      'Color','r');
        legend('Motif A','Motif B (unscaled)', ...
       'Location','best', ...
       'Box','on');

        subplot(2,1,2);
        plot(motifA,'b','LineWidth',2); hold on;
        plot(motifB_scaled,'g','LineWidth',2);
        title(sprintf('AFTER Scaling at %.3fx (Distance: %.2f)', ...
      bestScale, dist_after), 'Color',[0 0.6 0]);
        legend('Motif A','Motif B (scaled)', ...
       'Location','best', ...
       'Box','on');
    end

    % ----- return struct -----
    result = struct();
    result.scales     = scales;
    result.distances  = distances;
    result.bestScale  = bestScale;
    result.bestDist   = bestDistFinal;
    result.bestIdxA   = bestIdxA;
    result.bestIdxB_scaled   = bestIdxB_scaled;
    result.bestIdxB_unscaled = bestIdxB_unscaled;
    result.motifLength = m;
    result.minlag = minlag;
end

% ===================== helpers =====================

function y = rescale_linear(x, sf)
    n = length(x);
    new_len = max(2, round(n*sf));
    y = interp1(1:n, x, linspace(1,n,new_len), 'linear')';
end

function w = extract_window(x, startIdx, m)
    startIdx = max(1, round(startIdx));
    endIdx = min(startIdx+m-1, length(x));
    w = x(startIdx:endIdx);
    if length(w) < m
        w = [w; zeros(m-length(w),1)];
    end
end