%% Performances in function of the number of snapshots

clear variables
load damasdemo
close all

%% Parameters of the simulation

nbSources=4;

snr = -10;

SNAPS = [1 2 5 10 20 50 100 200 500 1000 2000 5000];

k = 80;
nbMethods = 8;

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

        
        % SWF-COMET1
        
        [XSP1, amps_SP1, ~, ~, amps_SP1_pinv] = sfw_comet1(Pmic, k, DD, XX, 0, 0, nbSources, [LBx LBy LBz]-0.01, [UBx UBy UBz]+0.01);
        SOURCES_POS_EST{ii, s, 2} = XSP1;
        SOURCES_AMPS_EST{ii, s, 2} = amps_SP1;
        % least-squares reestimation
        SOURCES_AMPS_EST{ii, s, 9} = amps_SP1_pinv;
        
        % SWF-COMET2       
        [XSP2, amps_SP2, ~, ~, amps_SP2_pinv] = sfw_comet2(Pmic, k, DD, XX, 0, 0, nbSources, [LBx LBy LBz]-0.01, [UBx UBy UBz]+0.01);
        SOURCES_POS_EST{ii, s, 3} = XSP2;
        SOURCES_AMPS_EST{ii, s, 3} = amps_SP2;
        
        SOURCES_AMPS_EST{ii, s, 10} = amps_SP2_pinv;
        
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
        
                
        % HR-CLEAN-SC
        [Xcsc, Pcsc, H] = hr_clean_sc_dr(DD,nbSources,XX, Pmic, k, 1, 0.25);
        SOURCES_POS_EST{ii, s, 8} = Xcsc;
        SOURCES_AMPS_EST{ii, s, 8} = Pcsc;
    end
end

save data_fig_SNAPS

%%

load data_fig_SNAPS

abss = SNAPS;
xtitle = 'Number of snapshots';

plot_perf(SOURCES_POS, SOURCES_AMPS, SOURCES_POS_EST, SOURCES_AMPS_EST, nbSources, abss, xtitle)

