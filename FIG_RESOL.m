%% Resolution

clear variables
load damasdemo
close all

%% Parameters of the simulation
nbSources=2;

% tested coordinate
coord = 'Z';

if coord == 'Z'
    % frequency
    k = 500;
    
    % index of the coordinate
    idxc = 3;
    
    % space between sources
    DELTAS = linspace(0.01, 0.2, 200);
    
else
    k = 50;
    
    idxc = 1;
    DELTAS = linspace(0.01, 0.2, 200);
end

snr = -10;

snapshots = 500;
nbMethods = 7;


nbPlots = length(DELTAS);

f = 340 * k /2/pi;

sample_size = 1;

SOURCES_POS = cell(nbPlots, sample_size);
SOURCES_AMPS = cell(nbPlots, sample_size);

SOURCES_POS_EST = cell(nbPlots, sample_size, nbMethods);
SOURCES_AMPS_EST = cell(nbPlots, sample_size, nbMethods);

DATAS = cell(nbPlots, sample_size);


Niter = 50;
% domaine

if coord == 'X'
    
    LBx = -0.5;
    UBx = 0.5;
    LBy = -0.02;
    UBy = 0.02;
    LBz = 3.98;
    UBz = 4.02;
    
else
    
    LBx = -0.02;
    UBx = 0.02;
    LBy = -0.02;
    UBy = 0.02;
    LBz = 3.5;
    UBz = 4.5;
    
end

step = 0.01;

% discretisation
xx = (LBx:step:UBx)';
yy = (LBy:step:UBy)';
zz = (LBz:step:UBz)';

[Xg, Yg, Zg] = meshgrid(xx, yy, zz);

XX = [Xg(:) Yg(:) Zg(:)];

% power
AA = 1;

SS = cell(nbMethods,1);

z = 0;


for ii = 1:nbPlots
    
    delta = DELTAS(ii);
    for s = 1:sample_size
        z = z+1;
        waitbar(z/(sample_size*nbPlots))
        
        if coord == 'Z'
            XS = [0 0 4-delta ; 0 0 4+delta];
        else
            XS = [-delta 0 4 ; delta 0 4];
        end
        
        a = randn(2,snapshots) + 1i * randn(2,snapshots);
        
        a = a .* sqrt(AA/2);
        AAgt = mean(abs(a).^2, 2);
        
        Dom = dictionary(Pmic, XS, k);
        Y0 = Dom * a;
        
        noise = randn(size(Y0));
        
        noise = noise / norm(noise, 'fro') * norm(Y0, 'fro') * 10^(-snr/20);
        Y = Y0 + noise;
        
        DD = Y * Y' / snapshots;
        
        DATAS{ii, s} = DD;
        
        SOURCES_POS{ii, s} = XS;
        SOURCES_AMPS{ii, s} = AAgt;
        
        
        %         % SFW-CMF
        %         [XSFWm, amps_SFWm] = sfw_cmf(Pmic, k, DD, XX, 0, 0, 1e-6, nbSources, [LBx LBy LBz]-0.01, [UBx UBy UBz]+0.01);
        %         SOURCES_POS_EST{ii, s, 1} = XSFWm;
        %         SOURCES_AMPS_EST{ii, s, 1} = amps_SFWm;
        %
        %         SS{1} = [SS{1} ; XSFWm(:, idxc) delta*ones(size(XSFWm, 1), 1) amps_SFWm];
        
        %         % SWF-SPICE
        %
        %         [XSP1, amps_SP1] = sfw_spice(Pmic, k, DD, XX, 0, 0, nbSources, [LBx LBy LBz]-0.01, [UBx UBy UBz]+0.01);
        %          SOURCES_POS_EST{ii, s, 2} = XSP1;
        %          SOURCES_AMPS_EST{ii, s, 2} = amps_SP1;
        %
        %         SS{2} = [SS{2} ; XSP1(:, idxc) delta*ones(size(XSP1, 1), 1) amps_SP1];
        
        % SWF-SPICE2
        
        [XSP2, amps_SP2] = sfw_spice2(Pmic, k, DD, XX, 0, 0, nbSources, [LBx LBy LBz]-0.01, [UBx UBy UBz]+0.01);
        SOURCES_POS_EST{ii, s, 3} = XSP2;
        SOURCES_AMPS_EST{ii, s, 3} = amps_SP2;
        
        SS{3} = [SS{3} ; XSP2(:, idxc) delta*ones(size(XSP2, 1), 1) amps_SP2];
        
        
        [XM, Pmest] = MUSIC_local_1D(DD, nbSources,XX,Pmic, k);
        SOURCES_POS_EST{ii, s, 4} = XM;
        SOURCES_AMPS_EST{ii, s, 4} = Pmest;
        
        SS{4} = [SS{4} ; XM(:, idxc) delta*ones(size(XM, 1), 1) Pmest];
        
        
        [xobf, q_OBF] = OBF(DD,nbSources,XX, Pmic, k);
        SOURCES_POS_EST{ii, s, 5} = xobf;
        SOURCES_AMPS_EST{ii, s, 5} = q_OBF;
        
        SS{5} = [SS{5} ; xobf(:, idxc) delta*ones(size(xobf, 1), 1) q_OBF];
        
        
        %         % OMP
        %          [xomp, q_OMP] = OMPDAMAS_DR(DD,nbSources,XX, Pmic, k);
        %          SOURCES_POS_EST{ii, s, 6} = xomp;
        %          SOURCES_AMPS_EST{ii, s, 6} = q_OMP;
        %
        %          SS{6} = [SS{6} ; xomp(:, idxc) delta*ones(size(xomp, 1), 1) q_OMP];
        
        
        
        % CLEAN-SC
        [Xest, P, H] = clean_sc_dr(DD,nbSources,XX, Pmic, k, 1);
        SOURCES_POS_EST{ii, s, 7} = Xest;
        SOURCES_AMPS_EST{ii, s, 7} = P;
        
        SS{7} = [SS{7} ; Xest(:, idxc) delta*ones(size(Xest, 1), 1) P'];
        
    end
end


save(['resol' coord])

%%
if coord == 'Z'
    load resolZ
else
    load resolX
end
C = colororder();

scatter(SS{3}(:, 2)*2, SS{3}(:, 1), SS{3}(:, 3)*50+eps, C(3, :), 'o')
hold on

scatter(SS{4}(:, 2)*2, SS{4}(:, 1), SS{4}(:, 3)*50+eps, C(4, :), 's')
scatter(SS{5}(:, 2)*2, SS{5}(:, 1), SS{5}(:, 3)*50+eps, C(5, :), '^')
scatter(SS{7}(:, 2)*2, SS{7}(:, 1), max(SS{7}(:, 3),0)*100+eps, C(7, :), '+')


legend('SP2', 'MUSIC', 'OBF', 'CLEAN-SC')

xlabel('\delta')
if coord == 'X'
    ylabel('Estimated X coordinate')
else
    ylabel('Estimated Z coordinate')
end


