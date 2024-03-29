%% SFW for source localization

% Demo with experimental measurements

clear all

% load the data
load data_sfw
close all

%% Parameters of the simulation

% do we test differential evolution (long)
DE = true;

% frequency
f = 340 * k /2/pi;

% domain of interest
LBx = -2;
UBx = 1;
LBy = -1;
UBy = 0;
LBz = 4;
UBz = 5;

% discretization step (initialization grid for the local optimizations in
% SFW)
step = 0.05;

% discretisation
xx = (LBx:step:UBx)';
yy = (LBy:step:UBy)';
zz = (LBz:step:UBz)';
[Xg, Yg, Zg] = meshgrid(xx, yy, zz);
XX = [Xg(:) Yg(:) Zg(:)];

% number of iterations for SFW
nbSources = 8;

% beamforming grid
stepb = 0.01;
xxb = (LBx:stepb:UBx)';
yyb = (LBy:stepb:UBy)';
zzb = 4.6;
[Xgb, Ygb, Zgb] = meshgrid(xxb, yyb, zzb);
XXb = [Xgb(:) Ygb(:) Zgb(:)];

% steering vectors
dic = dictionary(Pmic, XXb, k);
dic = dic ./ sqrt(sum(abs(dic.^2), 1));

%% noiseless case

% covariance matrix
DD = data*data' / size(data, 2);%Data;

% beamforming map
beam = real(sum(conj(dic) .* (DD * dic)));

% SFW-CMF
%[XSFWm, amps_SFWm] = sfw_cmf(Pmic, k, DD, XX, 0, 0, 1e-6, nbSources, [LBx LBy LBz]-0.01, [UBx UBy UBz]+0.01);

% SWF-COMET1
%[XSP1, amps_SP1] = sfw_comet1(Pmic, k, DD, XX, 0, 0, nbSources, [LBx LBy LBz]-0.01, [UBx UBy UBz]+0.01);

% SWF-COMET2
tic
[XSP2, amps_SP2, sn] = sfw_comet2(Pmic, k, DD, XX, 0, 0, nbSources, [LBx LBy LBz]-0.01, [UBx UBy UBz]+0.01);
TSP2 = toc;

% MUSIC
tic
[XM, Pmest] = MUSIC_local_exp(DD, nbSources,XX,Pmic, k);
TM = toc;
%  [xobf, q_OBF] = OBF(DD,nbSources,XX, Pmic, k);

% OMP
[Xomp, q_OMP] = OMPDAMAS_DR(DD,nbSources,XX, Pmic, k);

% CLEAN-SC
tic
[Xcsc, q_CSC] = clean_sc_dr(DD,nbSources,XX, Pmic, k, 1);
TCSC = toc;

tic
[Xhrcsc, q_HRCSC] = hr_clean_sc(DD,nbSources,XX, Pmic, k, 1, 0.25);
THRCSC = toc;

% differential evolution
if DE
    
    params = getdefaultparams();
    params.NP = 128*4;
    params.F = 0.8;
    params.CR = 0.55;
    params.maxiter = 200;
    params.displayResults = 0;
    variables = {
        'X1', [-2, 1 ; -1 0 ; 4 5], [0 ;0 ;0], (([LBx LBy LBz] + [UBx UBy UBz])/2 + (rand(1, 3) - 0.5)*2)' ;
        'X2', [-2, 1 ; -1 0 ; 4 5], [0 ;0 ;0], (([LBx LBy LBz] + [UBx UBy UBz])/2 + (rand(1, 3) - 0.5)*2)' ;
        'X3', [-2, 1 ; -1 0 ; 4 5], [0 ;0 ;0], (([LBx LBy LBz] + [UBx UBy UBz])/2 + (rand(1, 3) - 0.5)*2)' ;
        'X4', [-2, 1 ; -1 0 ; 4 5], [0 ;0 ;0], (([LBx LBy LBz] + [UBx UBy UBz])/2 + (rand(1, 3) - 0.5)*2)' ;
        'p1', [0, 1e6], 0, 3 ;
        'p2', [0, 1e6], 0, 3 ;
        'p3', [0, 1e6], 0, 3 ;
        'p4', [0, 1e6], 0, 3 ;
        'noise', [0, 1e6], 0, 1
        }
    
    settings.Pmic = Pmic;
    settings.k = k;
    settings.DD = DD;
    tic
    [bestmem, bestval] = differentialevolution(params, variables, @objde, settings);
    
    depos = reshape(bestmem(1:12), 3, 4);
    depow = bestmem(13:16);
    
    TDE = toc
end
save exp1

%% noiseless case - plots
load exp1

% ground truth
load Xgt

C = colororder();

% beamforming
imagesc(xxb, yyb, 10*log10(reshape(beam, 101, 301)))
axis xy
axis image
colormap(hot)
hold on
scatter(Xgt(:, 1), Xgt(:, 2), 100, 'o', 'filled', 'MarkerFaceColor', [0.8 0.8 0.8])
xlabel('X (m)')
ylabel('Y (m)')
xlim([-2, 1])
ylim([-1, 0])
cb = colorbar;
ylabel(cb,'Power (dB)')%,'FontSize',16,'Rotation',270);

% COMET2 and MUSIC

figure('Position', [100, 100, 1200, 500])

% XY
subplot(2, 3, 1)
scatter(Xgt(:, 1), Xgt(:, 2), 200, 'o', 'filled', 'MarkerFaceColor', [0.8 0.8 0.8])
hold on
scatter(XSP2(:, 1),XSP2(:, 2), amps_SP2/1e3+eps, C(3, :), '+', 'linewidth', 3)
scatter(XM(:, 1),XM(:, 2),Pmest/1e3+eps, C(7, :), 'x', 'linewidth', 3)
axis image
legend('Sources', 'COMET2', 'MUSIC')
xlabel('X (m)')
ylabel('Y (m)')
xlim([-2, 1])
ylim([-1, 0])

% XZ
subplot(2, 3, 4)
scatter(Xgt(:, 1), Xgt(:, 3), 200, 'o', 'filled', 'MarkerFaceColor', [0.8 0.8 0.8])
hold on
scatter(XSP2(:, 1),XSP2(:, 3), amps_SP2/1e3+eps, C(3, :), '+', 'linewidth', 3)
scatter(XM(:, 1), XM(:, 3),Pmest/1e3+eps, C(7, :), 'x', 'linewidth', 3)
axis image
xlabel('X (m)')
ylabel('Z (m)')
xlim([-2, 1])
ylim([4, 5])

% CLEAN-SC and differential evolution

% XY
subplot(2, 3, 2)
scatter(Xgt(:, 1), Xgt(:, 2), 200, 'o', 'filled', 'MarkerFaceColor', [0.8 0.8 0.8])
hold on
scatter(Xcsc(:, 1),Xcsc(:, 2),q_CSC/1e3+eps, C(2, :), '+', 'linewidth', 3)
scatter(Xhrcsc(:, 1),Xhrcsc(:, 2),q_HRCSC/1e3+eps, C(6, :), 'x', 'linewidth', 3)
    legend('Sources','CLEAN-SC', 'HR-CLEAN-SC') 
axis image

xlabel('X (m)')
ylabel('Y (m)')
xlim([-2, 1])
ylim([-1, 0])

% YZ


 if DE
subplot(2, 3, 5)
scatter(Xgt(:, 1), Xgt(:, 3), 200, 'o', 'filled', 'MarkerFaceColor', [0.8 0.8 0.8])

hold on
scatter(Xcsc(:, 1),Xcsc(:, 3),q_CSC/1e3+eps, C(2, :), '+', 'linewidth', 3)
scatter(Xhrcsc(:, 1),Xhrcsc(:, 3),q_HRCSC/1e3+eps, C(6, :), 'x', 'linewidth', 3)
axis image

xlabel('X (m)')
ylabel('Y (m)')
xlim([-2, 1])
ylim([4, 5])
% XY
subplot(2, 3, 3)
scatter(Xgt(:, 1), Xgt(:, 2), 200, 'o', 'filled', 'MarkerFaceColor', [0.8 0.8 0.8])
hold on

    scatter(depos(1, :),depos(2, :),depow/1e3+eps, C(5, :), 's', 'linewidth', 3)
    legend('Sources', 'DE')

axis image

xlabel('X (m)')
ylabel('Y (m)')
xlim([-2, 1])
ylim([-1, 0])

% YZ
subplot(2, 3, 6)
scatter(Xgt(:, 1), Xgt(:, 3), 200, 'o', 'filled', 'MarkerFaceColor', [0.8 0.8 0.8])

hold on

    scatter(depos(1, :),depos(3, :),depow/1e3+eps, C(5, :), 's', 'linewidth', 3)

axis image
xlabel('X (m)')
ylabel('Z (m)')
xlim([-2, 1])
ylim([4, 5])
 end
%% noisy case

% noise
SNR = -10;

Pdata = var(data(:));
Pnoise = Pdata * 10^(-SNR/10);

noise = (randn(size(data)) + 1i * randn(size(data))) * sqrt(Pnoise)/sqrt(2);
datan = data + noise;

% covariance matrix
DD = datan*datan' / size(data, 2);



% % SFW-CMF
%  [XSFWm, amps_SFWm] = sfw_cmf(Pmic, k, DD, XX, 0, 0, 1e-6, nbSources, [LBx LBy LBz]-0.01, [UBx UBy UBz]+0.01);
%
%  % SWF-COMET1
%
%  [XSP1, amps_SP1] = sfw_comet1(Pmic, k, DD, XX, 0, 0, nbSources, [LBx LBy LBz]-0.01, [UBx UBy UBz]+0.01);


% SWF-COMET2
[XSP2, amps_SP2] = sfw_comet2(Pmic, k, DD, XX, 0, 0, nbSources, [LBx LBy LBz]-0.01, [UBx UBy UBz]+0.01);

% MUSIC
[XM, Pmest] = MUSIC_local_exp(DD, nbSources,XX,Pmic, k);

%[xobf, q_OBF] = OBF(DD,nbSources,XX, Pmic, k);

% OMP
%[xomp, q_OMP] = OMPDAMAS_DR(DD,nbSources,XX, Pmic, k);

% CLEAN-SC
[Xcsc, q_CSC] = clean_sc_dr(DD,nbSources,XX, Pmic, k, 1);

[Xhrcsc, q_HRCSC] = hr_clean_sc_dr(DD,nbSources,XX, Pmic, k, 1, 0.25);

% DE
if DE
    
    params = getdefaultparams();
    params.NP = 128*4;
    params.F = 0.8;
    params.CR = 0.55;
    params.maxiter = 200;
    params.displayResults = 0;
    variables = {
        'X1', [-2, 1 ; -1 0 ; 4 5], [0 ;0 ;0], (([LBx LBy LBz] + [UBx UBy UBz])/2 + (rand(1, 3) - 0.5)*2)' ;
        'X2', [-2, 1 ; -1 0 ; 4 5], [0 ;0 ;0], (([LBx LBy LBz] + [UBx UBy UBz])/2 + (rand(1, 3) - 0.5)*2)' ;
        'X3', [-2, 1 ; -1 0 ; 4 5], [0 ;0 ;0], (([LBx LBy LBz] + [UBx UBy UBz])/2 + (rand(1, 3) - 0.5)*2)' ;
        'X4', [-2, 1 ; -1 0 ; 4 5], [0 ;0 ;0], (([LBx LBy LBz] + [UBx UBy UBz])/2 + (rand(1, 3) - 0.5)*2)' ;
        'p1', [0, 1e6], 0, 3 ;
        'p2', [0, 1e6], 0, 3 ;
        'p3', [0, 1e6], 0, 3 ;
        'p4', [0, 1e6], 0, 3 ;
        'noise', [0, 1e6], 0, 1
        }
    
    settings.Pmic = Pmic;
    settings.k = k;
    settings.DD = DD;
    tic
    [bestmem, bestval] = differentialevolution(params, variables, @objde, settings);
    
    depos = reshape(bestmem(1:12), 3, 4);
    depow = bestmem(13:16);
    
    TDE = toc
end
save exp2

%% noisy case - plots
load exp2
load Xgt

figure('Position', [100, 100, 1200, 500])

subplot(2, 3, 1)
scatter(Xgt(:, 1), Xgt(:, 2), 100, 'o', 'filled', 'MarkerFaceColor', [0.8 0.8 0.8])
hold on
scatter(XSP2(:, 1),XSP2(:, 2), amps_SP2/1e3+eps,C(3, :), '+', 'linewidth', 3)
scatter(XM(:, 1),XM(:, 2),Pmest/1e3+eps, C(7, :),'x', 'linewidth', 3)
legend('Sources', 'COMET2', 'MUSIC')
xlabel('X (m)')
ylabel('Y (m)')
axis image
xlim([-2, 1])
ylim([-1, 0])

subplot(2, 3, 4)
scatter(Xgt(:, 1), Xgt(:, 3), 100, 'o', 'filled', 'MarkerFaceColor', [0.8 0.8 0.8])
hold on
scatter(XSP2(:, 1),XSP2(:, 3), amps_SP2/1e3+eps,C(3, :), '+', 'linewidth', 3)
scatter(XM(:, 1),XM(:, 3),Pmest/1e3+eps, C(7, :),'x', 'linewidth', 3)
xlabel('X (m)')
ylabel('Z (m)')
axis image
xlim([-2, 1])
ylim([4, 5])

subplot(2, 3, 2)
scatter(Xgt(:, 1), Xgt(:, 2), 100, 'o', 'filled', 'MarkerFaceColor', [0.8 0.8 0.8])
hold on
scatter(Xcsc(:, 1),Xcsc(:, 2),max(q_CSC,0)/1e3+eps, C(2, :), '+', 'linewidth', 3)
scatter(Xhrcsc(:, 1),Xhrcsc(:, 2),max(q_HRCSC, 0)/1e3+eps, C(6, :), 'x', 'linewidth', 3)


    legend('Sources', 'CLEAN-SC', 'HR-CLEAN-SC')


xlabel('X (m)')
ylabel('Y (m)')
axis image
xlim([-2, 1])
ylim([-1, 0])

subplot(2, 3, 5)
scatter(Xgt(:, 1), Xgt(:, 3), 100, 'o', 'filled', 'MarkerFaceColor', [0.8 0.8 0.8])
hold on
scatter(Xcsc(:, 1),Xcsc(:, 3),max(q_CSC,0)/1e3+eps,C(2, :),  '+', 'linewidth', 3)
scatter(Xhrcsc(:, 1),Xhrcsc(:, 3),max(q_HRCSC, 0)/1e3+eps, C(6, :), 'x', 'linewidth', 3)



xlabel('X (m)')
ylabel('Z (m)')
axis image
xlim([-2, 1])
ylim([4, 5])

if DE
    
subplot(2, 3, 3)
scatter(Xgt(:, 1), Xgt(:, 2), 100, 'o', 'filled', 'MarkerFaceColor', [0.8 0.8 0.8])
hold on

    scatter(depos(1, :),depos(2, :),depow/1e3+eps, C(5, :), 's', 'linewidth', 3)
    legend('Sources', 'DE')    


xlabel('X (m)')
ylabel('Y (m)')
axis image
xlim([-2, 1])
ylim([-1, 0])

subplot(2, 3, 6)
scatter(Xgt(:, 1), Xgt(:, 3), 100, 'o', 'filled', 'MarkerFaceColor', [0.8 0.8 0.8])
hold on

    scatter(depos(1, :),depos(3, :),depow/1e3+eps, C(5, :), 's', 'linewidth', 3)


xlabel('X (m)')
ylabel('Z (m)')
axis image
xlim([-2, 1])
ylim([4, 5])

end
