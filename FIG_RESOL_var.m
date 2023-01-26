%% Influcence of the noise on precision

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
    DELTAS = 0.01:0.01:0.07;
else
    k = 50;

    idxc = 1;
    DELTAS = 0.01:0.01:0.07;
    
end

snr = -10;

snapshots = 500;
nbMethods = 2;


nbPlots = length(DELTAS);


f = 340 * k /2/pi;

sample_size = 100;

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


        % SWF-COMET2

        [XSP2, amps_SP2] = sfw_comet2(Pmic, k, DD, XX, 0, 0, nbSources, [LBx LBy LBz]-0.01, [UBx UBy UBz]+0.01);
         SOURCES_POS_EST{ii, s, 1} = XSP2;
         SOURCES_AMPS_EST{ii, s, 1} = amps_SP2;
                                  
      
         [Xest, P, H] = hr_clean_sc_dr(DD,nbSources,XX, Pmic, k, 1, 0);
         SOURCES_POS_EST{ii, s, 2} = Xest;
         SOURCES_AMPS_EST{ii, s, 2} = P;   
         
    end
end


save(['resol_var' coord])

%%
if coord == 'Z'
load resol_varZ
Z = 4;
else
    load resol_varX
    Z = 0;
end



 XSP2 = zeros(nbPlots, 2, sample_size);
 XCSC = zeros(nbPlots, 2, sample_size);

 ESP2 = zeros(nbPlots, 2, sample_size);
 ECSC = zeros(nbPlots, 2, sample_size);



for ii = 1:nbPlots


    for s = 1:sample_size
        XSP2(ii, :, s) = sort(SOURCES_POS_EST{ii, s, 1}(:, idxc));
        XCSC(ii, :, s) = sort(SOURCES_POS_EST{ii, s, 2}(:, idxc));
        
        ESP2(ii, 1, s) = XSP2(ii, 1, s) + DELTAS(ii) - Z;
        ESP2(ii, 2, s) = XSP2(ii, 2, s) - DELTAS(ii) - Z;
        ECSC(ii, 1, s) = XCSC(ii, 1, s) + DELTAS(ii) - Z;
        ECSC(ii, 2, s) = XCSC(ii, 2, s) - DELTAS(ii) - Z;

    end
end
    BSP2 = mean(ESP2(:, :, :), 3);
    BCSC = mean(ECSC(:, :, :), 3);


    VSP2 = var(XSP2(:, :, :),[], 3);
    VCSC = var(XCSC(:, :, :),[], 3);


    

C = colororder();

LW = 2;

MSP2 = mean(ESP2(:, :, :), 3);
Q1SP2 = quantile(ESP2, 0.25, 3);
Q3SP2 = quantile(ESP2, 0.75, 3);

MCSC = mean(ECSC(:, :, :), 3);
Q1CSC = quantile(ECSC, 0.25, 3);
Q3CSC = quantile(ECSC, 0.75, 3);

figure('Position', [100, 100, 800, 600])


subplot(2, 1, 1)
errorbar(2*DELTAS-0.001, MSP2(:, 1), -(Q1SP2(:, 1) - MSP2(:, 1)),  (Q3SP2(:, 1) - MSP2(:, 1)), '+', 'Color', C(3, :), 'linewidth', LW)
hold on
errorbar(2*DELTAS+0.001, MCSC(:, 1), -(Q1CSC(:, 1) - MCSC(:, 1)),  (Q3CSC(:, 1) - MCSC(:, 1)), '^', 'Color', C(6, :), 'linewidth', LW)

plot(DELTAS*4-0.04, zeros(size(DELTAS)), '-k')
xlim(2*[0.015, 0.075])
ylim([-0.05, 0.02])
legend('COMET2', 'HRCSC')
if coord == 'Z'
    ylabel('Error (m) (Z < 4m)')
else
    ylabel('Error (m) (X < 0m)')
end
xlabel('\delta')

subplot(2, 1, 2)
errorbar(2*DELTAS-0.001, MSP2(:, 2), -(Q1SP2(:, 2) - MSP2(:, 2)),  (Q3SP2(:, 2) - MSP2(:, 2)), '+', 'Color', C(3, :), 'linewidth', LW)
hold on
errorbar(2*DELTAS+0.001, MCSC(:, 2), -(Q1CSC(:, 2) - MCSC(:, 2)),  (Q3CSC(:, 2) - MCSC(:, 2)), '^', 'Color', C(6, :), 'linewidth', LW)

plot(DELTAS*4-0.04, zeros(size(DELTAS)), '-k')
  xlim(2*[0.015, 0.075])
ylim([-0.02, 0.04]) 
if coord == 'Z'
    ylabel('Error (m) (Z > 4m)')
else
    ylabel('Error (m) (X > 0m)')
end
xlabel('\delta')




