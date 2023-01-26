function plot_perf(SOURCES_POS, SOURCES_AMPS, SOURCES_POS_EST, SOURCES_AMPS_EST, nbSources, abss, xtitle)



if strcmp(xtitle, 'Number of snapshots')
    plotlog = @loglog;
else
    plotlog = @semilogy;
end

[nbPlots, sample_size, nbMethods] = size(SOURCES_POS_EST);

% errors
EPOS = zeros(nbPlots, nbMethods, sample_size);
EAMPS = zeros(nbPlots, nbMethods+2, sample_size);

% true is MUSIC has found all the sources
music_has_found_all = zeros(nbPlots, sample_size);

for ii = 1:nbPlots

    for s = 1:sample_size
        for m = 1:nbMethods

            [epos, eamps] = compute_errors(SOURCES_POS_EST{ii, s, m}, SOURCES_POS{ii, s}, SOURCES_AMPS_EST{ii, s, m}, SOURCES_AMPS{ii, s});
            EPOS(ii, m, s) = epos;
            EAMPS(ii, m, s) = eamps;
            
            if m == 4                
                music_has_found_all(ii, s) = (size(SOURCES_POS_EST{ii, s, m}, 1)  == nbSources);
            end
            
            
        end

        
        % power errors for the reestimated COMET
        m = 9;
        [epos, eamps] = compute_errors(SOURCES_POS_EST{ii, s, 2}, SOURCES_POS{ii, s}, SOURCES_AMPS_EST{ii, s, m}, SOURCES_AMPS{ii, s});
        EAMPS(ii, m, s) = eamps;
        m = 10;
        [epos, eamps] = compute_errors(SOURCES_POS_EST{ii, s, 3}, SOURCES_POS{ii, s}, SOURCES_AMPS_EST{ii, s, m}, SOURCES_AMPS{ii, s});
        EAMPS(ii, m, s) = eamps;
    end
end

% minimal ratio of identified sources for MUSIC
ratiom = sum(music_has_found_all, 2)/sample_size;
ratiomin = 0.9;

% linewidth
LW = 2;
% marker size
MS = 15;

% MSE
mEPOS = mean(EPOS, 3);
mEAMPS = mean(EAMPS, 3);

C = colororder();



figure('Position', [100, 100, 800, 600])

% position

subplot(2, 2, 1)
plotlog(abss,  mEPOS(:, 1), '-s', 'Color', C(1, :), 'linewidth', LW, 'markersize', MS)
hold on

if strcmp(xtitle, 'Number of snapshots')
    plot(abss, mEPOS(:, 2).* (abss >= 128)', '-x', 'Color', C(2, :), 'linewidth', LW, 'markersize', MS)
else
    plot(abss, mEPOS(:, 2), '-x', 'Color', C(2, :), 'linewidth', LW, 'markersize', MS)
end
plot(abss, mEPOS(:, 3), '-+', 'Color', C(3, :), 'linewidth', LW, 'markersize', MS)
plot(abss, mEPOS(:, 5), '--', 'Color', C(4, :), 'linewidth', LW, 'markersize', MS)

legend('CMF', 'COMET1', 'COMET2', 'OBF')

xlabel(xtitle)
ylabel('Position MSE (m^2)')

subplot(2, 2, 2)
plotlog(abss, mEPOS(:, 3), '-+', 'Color', C(3, :), 'linewidth', LW, 'markersize', MS)
hold on
plot(abss, mEPOS(:, 4).*(ratiom >= ratiomin), '-x', 'Color', C(7, :), 'linewidth', LW, 'markersize', MS)
plot(abss, mEPOS(:, 6), '-', 'Color', C(1, :), 'linewidth', LW, 'markersize', MS)
plot(abss, mEPOS(:, 7), '--', 'Color', C(2, :), 'linewidth', LW, 'markersize', MS)
plot(abss, mEPOS(:, 8), '-^', 'Color', C(6, :), 'linewidth', LW, 'markersize', MS)

legend('COMET2', 'MUSIC', 'OMP', 'CSC', 'HRCSC')

xlabel(xtitle)
ylabel('Position MSE (m^2)')


% power

subplot(2, 2, 3)
plotlog(abss,  mEAMPS(:, 1), '-s', 'Color', C(1, :), 'linewidth', LW, 'markersize', MS)
hold on
if strcmp(xtitle, 'Number of snapshots')
    plot(abss, mEAMPS(:, 2).* (abss >= 128)', '-x', 'Color', C(2, :), 'linewidth', LW, 'markersize', MS)
else
    plot(abss, mEAMPS(:, 2), '-x', 'Color', C(2, :), 'linewidth', LW, 'markersize', MS)
end
plot(abss, mEAMPS(:, 3), '-+', 'Color', C(3, :), 'linewidth', LW, 'markersize', MS)
plot(abss, mEPOS(:, 5), '--', 'Color', C(4, :), 'linewidth', LW, 'markersize', MS)
plot(abss, mEAMPS(:, 10), '--o', 'Color', C(5, :), 'linewidth', LW, 'markersize', MS)


legend('CMF', 'COMET1', 'COMET2', 'OBF', 'COMET2LS')

xlabel(xtitle)
ylabel('Power MSE (Pa^4)')

subplot(2, 2, 4)

plotlog(abss, mEAMPS(:, 4).*(ratiom >= ratiomin), '-x', 'Color', C(7, :), 'linewidth', LW, 'markersize', MS)
hold on
plot(abss, mEAMPS(:, 6), '-', 'Color', C(1, :), 'linewidth', LW, 'markersize', MS)
plot(abss, mEAMPS(:, 7), '--', 'Color', C(2, :), 'linewidth', LW, 'markersize', MS)
plot(abss, mEAMPS(:, 8), '-^', 'Color', C(6, :), 'linewidth', LW, 'markersize', MS)
plot(abss, mEAMPS(:, 10), '--o', 'Color', C(5, :), 'linewidth', LW, 'markersize', MS)

legend('MUSIC', 'OMP', 'CSC', 'HRCSC', 'COMET2LS')

xlabel(xtitle)
ylabel('Power MSE (Pa^4)')

