% This script kills fascists and plots timeseries of all flowband model
% output in a given directory at Diamond Hill, Lake Wellman, Magnis Valley,
% and Dubris Valley
addpath '/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton/Manuscript/Figures/scripts/'
LW_ind = 13;
MV_ind = 31;
DV_ind = 46;
load 'DH_DATA/Geochronology data/cosmo_data.mat'
runs_dir = ['/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton/Modeling/Koutnik',...
    ' model/DH_code/TEST_DH/output/sealev_ens/deformation_and_sliding/'];

runs_dir_content = struct2table(dir(runs_dir));
ens_runs = runs_dir_content.name(3:end);

load_constraints_for_scoring;
%initialize tables that will hold scores
RMS_scores = table;
RMS_scores.name = ens_runs;
RMS_scores.tau = 1000.*~cellfun('isempty',strfind(RMS_scores.name, '1000'))...
    + 2000.*~cellfun('isempty',strfind(RMS_scores.name, '2000'));
RMS_scores.ocfac = 0.5.*~cellfun('isempty',strfind(RMS_scores.name, 'ocfac0_5'))...
    + 1.*~cellfun('isempty',strfind(RMS_scores.name, 'ocfac1'));
RMS_scores.crhshelf = 1e-5.*~cellfun('isempty',strfind(RMS_scores.name, 'shelf5'))...
    + 1e-6.*~cellfun('isempty',strfind(RMS_scores.name, 'shelf6')) ...
    + 1e-7.*~cellfun('isempty',strfind(RMS_scores.name, 'shelf7'));
RMS_scores.LR04_sealev = 1-cellfun('isempty',strfind(RMS_scores.name, 'FORCEPLEIST'));

RMS_scores.DV = nan*ones(length(ens_runs),1);
RMS_scores.MV = nan*ones(length(ens_runs),1);
RMS_scores.LW = nan*ones(length(ens_runs), 1);
RMS_scores.DH = nan*ones(length(ens_runs), 1);
RMS_scores.total = nan*ones(length(ens_runs), 1);

for jj = 1:length(ens_runs)
        plotline = {'Color', [0.6 0.6 0.6]};
  
    if strfind(ens_runs{jj}, 'smooth9ka')
        plotline = {'k'};
    end
    
    load([runs_dir,'/', ens_runs{jj}], 't_P', 't_P2', ...
        'x_P', 'x_P2', 'S_P', 'S_P2', 'S_at_GL');
   
    figure(1)
    hold on
    DV_ax = gca;
    plot(-t_P2, S_P2(:, DV_ind), plotline{:}, 'linewidth', 1.25);
    
    figure(2)
    hold on
    MV_ax = gca;
    plot(-t_P2, S_P2(:, MV_ind), plotline{:}, 'linewidth', 1.25);
    
    figure(3)
    hold on
    LW_ax = gca;
    plot(-t_P2, S_P2(:, LW_ind), plotline{:}, 'linewidth', 1.25);

    figure(4)
    hold on
    DH_ax = gca;
    plot(-t_P2, S_at_GL-S_at_GL(end), plotline{:}, 'linewidth', 1.25);
   
%% Now score this model run against geochronology

model_DV = interp1(-t_P2, S_P2(:, DV_ind),DAN_UM_DV_BV_ages);
RMS_scores.DV(jj) = sqrt(sum((DAN_UM_DV_BV_elev - model_DV).^2)./length(DAN_UM_DV_BV_elev));

model_MV = interp1(-t_P2, S_P2(:, MV_ind), MV_ages);
RMS_scores.MV(jj) = sqrt(sum((MV_elev - model_MV).^2)./length(MV_elev));

model_LW = interp1(-t_P2, S_P2(:, LW_ind), LW_ages);
RMS_scores.LW(jj) = sqrt(sum((LW_elev - model_LW).^2)./length(LW_elev));

model_DH = interp1(-t_P, S_P(:, 1)-S_at_GL(end), DH_ages);
RMS_scores.DH(jj) = sqrt(sum((DH_elev - model_DH).^2)./length(DH_elev));

% RMS_scores.total(jj) = (RMS_scores.DAN(jj).*RMS_scores.DV(jj).*RMS_scores.MV(jj)...
%     .*RMS_scores.LW(jj).*RMS_scores.DH(jj)).^(1/5);

if strfind(ens_runs{jj}, 'smooth9ka')
    RMS_scores.LR04_sealev(jj) = NaN;
end

end
%% total model score following section 2.4.1 (approach [a]) from Pollard et al., 2016

DV50 = median(RMS_scores.DV(1:end-1));
MV50 = median(RMS_scores.MV(1:end-1));
LW50 = median(RMS_scores.LW(1:end-1));
DH50 = median(RMS_scores.DH(1:end-1));

RMS_scores.total = exp(-(RMS_scores.DV./DV50 +...
    RMS_scores.MV./MV50 + RMS_scores.LW./LW50 + RMS_scores.DH./DH50));


%% model score statistics for each parameter
tau1000_mean = mean(RMS_scores.total(RMS_scores.tau==1000));
tau1000_std = std(RMS_scores.total(RMS_scores.tau==1000));
tau2000_mean = mean(RMS_scores.total(RMS_scores.tau==2000));
tau2000_std = std(RMS_scores.total(RMS_scores.tau==2000));

ocfac0_5_mean = mean(RMS_scores.total(RMS_scores.ocfac==0.5));
ocfac0_5_std = std(RMS_scores.total(RMS_scores.ocfac==0.5));
ocfac1_mean = mean(RMS_scores.total(RMS_scores.ocfac==1));
ocfac1_std = std(RMS_scores.total(RMS_scores.ocfac==1));

shelf5_mean = mean(RMS_scores.total(RMS_scores.crhshelf == 1e-5));
shelf5_std = std(RMS_scores.total(RMS_scores.crhshelf == 1e-5));
shelf6_mean = mean(RMS_scores.total(RMS_scores.crhshelf == 1e-6));
shelf6_std = std(RMS_scores.total(RMS_scores.crhshelf == 1e-6));
shelf7_mean = mean(RMS_scores.total(RMS_scores.crhshelf == 1e-7));
shelf7_std = std(RMS_scores.total(RMS_scores.crhshelf == 1e-7));

FORCEPLEIST_mean = mean(RMS_scores.total(RMS_scores.LR04_sealev == true));
FORCEPLEIST_std = std(RMS_scores.total(RMS_scores.LR04_sealev == true));
SPRATT_mean = mean(RMS_scores.total(RMS_scores.LR04_sealev == false));
SPRATT_std = std(RMS_scores.total(RMS_scores.LR04_sealev == false));
%% plot extracted datapoints

figure(1)
    hold on
    title({'\it Dubris and Bibra Valleys,'; 'Danum Platform, Updog Mtn '}, 'FontWeight', 'normal')
    plot(DAN_UM_DV_BV_ages, DAN_UM_DV_BV_elev, 'ko', 'MarkerFaceColor', [0.8, 0.8, 0.8]);
%     plot(DV_BV_ages, DV_BV_elev, 'ko', 'MarkerFaceColor', [0.8, 0.8, 0.8]);
    set(gca, 'Fontsize', 15, 'Xticklabel', [0 5 10 15 20])
    set(gcf, 'Position', [202   438   292   260], 'PaperPositionMode', 'auto')
    ylabel('Elevation (m)')
    xlabel ('Age (kyr BP)')
    grid on
    
figure(2)
    hold on
    title('\it Magnis Valley ', 'FontWeight', 'normal')
    plot(MV_ages, MV_elev, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8]);
    set(gca, 'Fontsize', 15, 'Xticklabel', [0 5 10 15 20])
    set(gcf, 'Position', [202   438   292   260], 'PaperPositionMode', 'auto')
    ylabel('Elevation (m)')
    xlabel ('Age (kyr BP)')
    grid on
    
figure(3)
    hold on
    title('\it Lake Wellman ', 'FontWeight', 'normal')
    plot(LW_ages, LW_elev, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8]);
    set(gca, 'Fontsize', 15, 'Xticklabel', [0 5 10 15 20])
    set(gcf, 'Position', [202   438   292   260], 'PaperPositionMode', 'auto')
    ylabel('Elevation (m)')
    xlabel ('Age (kyr BP)')
    grid on
    
figure(4)
    hold on
    title('\it Diamond Hill ', 'FontWeight', 'normal')
    plot(DH_ages, DH_elev, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8])
    set(gca, 'Fontsize', 15, 'Xticklabel', [0 5 10 15 20])
    set(gcf, 'Position', [202   438   292   260], 'PaperPositionMode', 'auto')
    ylabel('Height above modern glacier (m)')
    xlabel ('Age (kyr BP)')
    grid on
    
%% get the five highest scores and plot them in a better color
highest_scores = sort(RMS_scores.total);
highest_scores = highest_scores(end-6:end-1);
colors = [118,42,131
175,141,195
231,212,232
217,240,211
127,191,123
27,120,55]./255;

for jj = 1:length(highest_scores);
highest_scores_index(jj) = find(RMS_scores.total == highest_scores(jj));

load([runs_dir,'/', ens_runs{highest_scores_index(jj)}], 't_P', 't_P2', ...
        'x_P', 'x_P2', 'S_P', 'S_P2', 'S_at_GL');
   

    
    figure(1)

    hold on
    DV_ax = gca;
    plot(-t_P2, S_P2(:, DV_ind), 'Color', colors(jj,:), 'linewidth', 2);
    
    figure(2)
    hold on
    plot(-t_P2, S_P2(:, MV_ind), 'Color', colors(jj,:),  'linewidth', 2);
    
    figure(3)
    hold on
    LW_ax = gca;
    plot(-t_P2, S_P2(:, LW_ind), 'Color', colors(jj,:), 'linewidth', 2);

    
    figure(4)
    hold on
    DH_ax = gca;
    plot(-t_P2, S_at_GL-S_at_GL(end), 'Color', colors(jj,:),  'linewidth', 2);

end

%% plot model score stats
figure(6)
subplot(1,4,1)
FaceColor = [188,189,220]./255;
bar([1, 0], [mean(RMS_scores.total(RMS_scores.LR04_sealev == 1)),...
    mean(RMS_scores.total(RMS_scores.LR04_sealev == 0))], 'FaceColor', FaceColor)
set(gca, 'XTickLabel', {'Spratt', 'LR04'}, 'Fontsize', 15)
xlabel('Sea-level forcing')
ylim([0 0.035])
xlim([-0.5 1.5])
ylabel('Mean model score')
grid on

subplot(1,4,2)
bar([1000, 2000], [mean(RMS_scores.total(RMS_scores.tau==1000)),...
    mean(RMS_scores.total(RMS_scores.tau==2000))], 'FaceColor', FaceColor)
set(gca, 'fontsize', 15, 'YTickLabel', {'' '' '' ''})
xlabel({'\tau: rebound'; 'response time (yr)'})
ylim([0 0.035])
xlim([500 2500])
grid on

subplot(1,4,3)
bar([0.5 1], [mean(RMS_scores.total(RMS_scores.ocfac==0.5))...
    mean(RMS_scores.total(RMS_scores.ocfac==1))], 'FaceColor', FaceColor); 
xlabel({'Sub-shelf'; 'melting factor'})
set(gca, 'fontsize', 15, 'Xtick', [0.5 1], 'YTickLabel', {'' '' '' ''})
ylim([0 0.035])
xlim([0.25 1.25])
grid on

subplot(1,4,4)
 bar([-5 -6 -7], [mean(RMS_scores.total(RMS_scores.crhshelf == 1e-5))...
     mean(RMS_scores.total(RMS_scores.crhshelf == 1e-6))...
     mean(RMS_scores.total(RMS_scores.crhshelf == 1e-7))], 'FaceColor', FaceColor)
set(gca, 'XTickLabel', {'10^{-7}', '10^{-6}', '10^{-5}'}, 'YTickLabel', {'' '' '' ''})
xlabel({'Basal sliding parameter'; 'beneath modern ice shelf'})
set(gca, 'fontsize', 15)
ylim([0 0.035])
xlim([-7.5 -4.5])
grid on

set(gcf, 'Position', [419   500   743   205], 'PaperPositionMode', 'auto')

%% now plot up cosmo and algae data. Naw, keep this commented

% location_strings={'DV', 'BV', 'DAN', 'UM', 'DH', 'MV_floor', 'MV_walls', 'LW'};
% grouping_strings={{'DV', 'BV'}, {'DAN', 'UM'},{'MV_floor'},{'MV_walls'}, {'LW'}, {'DH'}};  
% scaling_scheme = 'St';
% sample_type = 'erratics';
% plot_cosmo(data, location_strings, grouping_strings, scaling_scheme, sample_type);
% close 5

%% plot up highest-scoring runs in map view

highest_scoring_runs = regexprep(RMS_scores.name(highest_scores_index),...
    '.mat','','ignorecase'); %need to remove .mat from each one

addpath /Users/trevorhillebrand/Documents/Antarctica/PSUICE3D/Model_output/
addpath('/Users/trevorhillebrand/Documents/MATLAB/Toolboxes/AntarcticMappingTools_v5.00//AntarcticMappingTools');


for mm = 1:length(highest_scoring_runs)
    if strfind(highest_scoring_runs{mm}, 'SPRATT')
    [tmp] = nc2mat(['/Volumes/Chunk/PSUICE/NCfiles/Run125kyr30km_sealev/'...
        'SPRATT/Run25kyr20km_SPRATT/', highest_scoring_runs{mm}, '/fort.92.nc'], {'hs', 'maskwater'});
    elseif strfind(highest_scoring_runs{mm}, 'FORCEPLEIST')
    [tmp] = nc2mat(['/Volumes/Chunk/PSUICE/NCfiles/Run125kyr30km_sealev/'...
        'FORCEPLEIST/Run25kyr20km_FORCEPLEIST/', highest_scoring_runs{mm}, '/fort.92.nc'], {'hs', 'maskwater'});
 
    end
        
    eval([highest_scoring_runs{mm}, '= tmp;']);
    
    if mm == 1
        alatd = double(tmp.alatd);
        alond = double(tmp.alond);
        time = double(tmp.time);
    end

figure(7)


subplot(2,3,mm)
RossSea = antmap;
mapzoom('Ross Ice Shelf',  1600)
bedmap2('bed', 'colorbar', 'none');
cmap_max = [8 46 107]./255;
cmap_min = [0.8 0.8 0.8];

cmap = ones(100,3);
cmap(:,1) = linspace(cmap_max(1),1, 100);
cmap(:,2) = linspace(cmap_max(2),1, 100);
cmap(:,3) = linspace(cmap_max(3),1, 100);
cmap(end,:) = cmap_min;
cmap(1,:) = [1 1 1]; %don't waste ink!

caxis([-1300 10])

colormap(cmap)



for plottime = [-12e3 -9e3, -8e3, -7e3, -6e3,-5e3, time(end)]
    plottime_ind = find(time==plottime);
    contourm(alatd, alond, double(tmp.maskwater(:,:,plottime_ind)), 1, 'Color', colors(mm,:), 'linewidth', 2)
    
end
 contourm(alatd, alond, double(tmp.maskwater(:,:,end)), 1, 'k', 'linewidth', 2)
    
bedmap2 ('gl', 'k--', 'linewidth', 1.5)
end
subplot(2,3,5)

cbar = colorbar;
ylabel(cbar, 'Bed Elevation (m)')
set(cbar, 'fontsize', 14, 'orientation', 'horizontal', 'location', 'southoutside')

set(gcf, 'Position', [130    98   988   607], 'PaperPositionMode', 'auto')
