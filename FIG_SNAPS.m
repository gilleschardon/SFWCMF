%% Performances in function of the number of snapshots

clear variables
load damasdemo
close all

%% Parameters of the simulation

nbSources=4;

snr = -10;

SNAPS = [1 2 5 10 20 50 100 200 500 1000 2000 5000];

k = 80;
nbMethods = 7;

nbPlots = length(SNAPS);


f = 340 * k /2/pi;

sample_size = 100;

SOURCES_POS = cell(nbPlots, sample_size);
SOURCES_AMPS = cell(nbPlots, sample_size);

SOURCES_POS_EST = cell(nbPlots, sample_size, nbMethods);
SOURCES_AMPS_EST = cell(nbPlots, sample_size, nbMethods+2);

DATAS = cell(nbPlots, sample_size);

% domain
LBx = -1;
UBx = 1;
LBy = -1;
UBy = 1;
LBz = 3;
UBz = 5;



step = 0.05;

% discretisation
xx = (LBx:step:UBx)';
yy = (LBy:step:UBy)';
zz = (LBz:step:UBz)';

[Xg, Yg, Zg] = meshgrid(xx, yy, zz);

XX = [Xg(:) Yg(:) Zg(:)];

% powers
AA = [1 1 0.5 0.1]';
% positions
XS=[0.01 -0.13 4.12 ; 0.11 -0.03 3.922 ; 0.5231 -0.3233 4.71 ; -0.434 0.322 3.43];

% for the waitbar
z = 0;
for ii = 1:nbPlots
snapshots = SNAPS(ii);
for s = 1:sample_size
        z = z+1;
        waitbar(z/(sample_size*nbPlots))
        
        % amplitudes at each snapshot
        a = randn(nbSources,snapshots) + 1i * randn(nbSources,snapshots);
        a = a .* sqrt(AA/2);
        
        AAgt = mean(abs(a).^2, 2);

        % dictionary of sources
        Dom = dictionary(Pmic, XS, k);
        
        % noiseless measurements
        Y0 = Dom * a;
        
        % noise
        noise = randn(size(Y0));

        noise = noise / norm(noise, 'fro') * norm(Y0, 'fro') * 10^(-snr/20);
        Y = Y0 + noise;
        
        % SCM
        DD = Y * Y' / snapshots;
       
        DATAS{ii, s} = DD;
        
        SOURCES_POS{ii, s} = XS;
        SOURCES_AMPS{ii, s} = AAgt;
        
        
        % SFW-CMF
        [XSFWm, amps_SFWm] = sfw_cmf(Pmic, k, DD, XX, 0, 0, 1e-6, nbSources, [LBx LBy LBz]-0.01, [UBx UBy UBz]+0.01);
        SOURCES_POS_EST{ii, s, 1} = XSFWm;
        SOURCES_AMPS_EST{ii, s, 1} = amps_SFWm;

        
        % SWF-SPICE
        
        [XSP1, amps_SP1, ~, ~, amps_SP1_pinv] = sfw_spice(Pmic, k, DD, XX, 0, 0, nbSources, [LBx LBy LBz]-0.01, [UBx UBy UBz]+0.01);
        SOURCES_POS_EST{ii, s, 2} = XSP1;
        SOURCES_AMPS_EST{ii, s, 2} = amps_SP1;
        % least-squares reestimation
        SOURCES_AMPS_EST{ii, s, 8} = amps_SP1_pinv;
        
        % SWF-SPICE2       
        [XSP2, amps_SP2, ~, ~, amps_SP2_pinv] = sfw_spice2(Pmic, k, DD, XX, 0, 0, nbSources, [LBx LBy LBz]-0.01, [UBx UBy UBz]+0.01);
        SOURCES_POS_EST{ii, s, 3} = XSP2;
        SOURCES_AMPS_EST{ii, s, 3} = amps_SP2;
        
        SOURCES_AMPS_EST{ii, s, 9} = amps_SP2_pinv;
        
        % MUSIC
        [XM, Pmest] = MUSIC_local(DD, nbSources,XX,Pmic, k);
        SOURCES_POS_EST{ii, s, 4} = XM;
        SOURCES_AMPS_EST{ii, s, 4} = Pmest;
        
        %OBF
        [xobf, q_OBF] = OBF(DD,nbSources,XX, Pmic, k);
        SOURCES_POS_EST{ii, s, 5} = xobf;
        SOURCES_AMPS_EST{ii, s, 5} = q_OBF;
        
        % OMP
        [xomp, q_OMP] = OMPDAMAS_DR(DD,nbSources,XX, Pmic, k);
        SOURCES_POS_EST{ii, s, 6} = xomp;
        SOURCES_AMPS_EST{ii, s, 6} = q_OMP;
        
        % CLEAN-SC
        [Xcsc, Pcsc, H] = clean_sc_dr(DD,nbSources,XX, Pmic, k, 1);
        SOURCES_POS_EST{ii, s, 7} = Xcsc;
        SOURCES_AMPS_EST{ii, s, 7} = Pcsc;
    end
end

save data_fig_SNAPS

%%

load data_fig_SNAPS



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
        
        % power errors for the reestimated SPICE
        m = 8;
        [epos, eamps] = compute_errors(SOURCES_POS_EST{ii, s, 2}, SOURCES_POS{ii, s}, SOURCES_AMPS_EST{ii, s, m}, SOURCES_AMPS{ii, s});
        EAMPS(ii, m, s) = eamps;
        m = 9;
        [epos, eamps] = compute_errors(SOURCES_POS_EST{ii, s, 3}, SOURCES_POS{ii, s}, SOURCES_AMPS_EST{ii, s, m}, SOURCES_AMPS{ii, s});
        EAMPS(ii, m, s) = eamps;
    end
end
figure

% x-axis
abss = SNAPS;

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


subplot(2, 1, 1)
loglog(abss,  mEPOS(:, 1), '-+', 'linewidth', LW, 'markersize', MS)
hold on
plot(abss, mEPOS(:, 2).* (SNAPS >= size(Pmic, 1))', '-x', 'linewidth', LW, 'markersize', MS)
plot(abss, mEPOS(:, 3), '-o', 'linewidth', LW, 'markersize', MS)
plot(abss, mEPOS(:, 4).*(ratiom >= ratiomin), '--s', 'linewidth', LW, 'markersize', MS)

plot(abss, mEPOS(:, 5), '--', 'linewidth', LW, 'markersize', MS)
plot(abss, mEPOS(:, 6), '-.', 'linewidth', LW, 'markersize', MS)
plot(abss, mEPOS(:, 7), '-', 'linewidth', LW, 'markersize', MS)
legend('CMF-N', 'SPICE1', 'SPICE2', 'MUSIC', 'OBF', 'OMP', 'CSC')

xlabel('snaps')
ylabel('Position MSE')

subplot(2, 1, 2)
loglog(abss, mEAMPS(:, 1), '-+', 'linewidth', LW, 'markersize', MS)
hold on
semilogy(abss, mEAMPS(:, 2).* (SNAPS >= size(Pmic, 1))', '-x', 'linewidth', LW, 'markersize', MS)
semilogy(abss, mEAMPS(:, 3), '-o', 'linewidth', LW, 'markersize', MS)
semilogy(abss, mEAMPS(:, 4).*(ratiom >= ratiomin), '--s', 'linewidth', LW, 'markersize', MS)
plot(abss, mEAMPS(:, 9), '--o', 'linewidth', LW, 'markersize', MS)


xlabel('snaps')
ylabel('Power MSE')
legend('CMF-N', 'SPICE1', 'SPICE2', 'MUSIC', 'SP2LS')

