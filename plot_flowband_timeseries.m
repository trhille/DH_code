% This script kills fascists and plots timeseries of all flowband model
% output in a given directory at Diamond Hill, Lake Wellman, Magnis Valley,
% and Dubris Valley
close all

addpath '/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton/Manuscript/Figures/scripts/'
LW_x_P_ind = 13;
MV_x_P_ind = 31;
DV_x_P_ind = 46;

% load constraints so you can score the model run
load 'DH_DATA/Geochronology data/geochron_constraints.mat'
load 'DH_DATA/Geochronology data/cosmo_data.mat'
load 'DH_DATA/Geochronology data/algae_data.mat'
% get the names of all runs in the ensemble you are scoring
runs_dir = ['/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton/Modeling/Koutnik',...
    ' model/DH_code/TEST_DH/output/deglaciation_scenarios/'];

runs_dir_content = struct2table(dir([runs_dir, '*.mat']));
ens_runs = runs_dir_content.name;

load([runs_dir,'/', ens_runs{1}], 't_P2')

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
RMS_scores.calv = 0.7.*~cellfun('isempty', strfind(RMS_scores.name, 'calv0_7')) ...
    + 1.*cellfun('isempty', strfind(RMS_scores.name, 'calv0_7'));
RMS_scores.LR04_sealev = 1-cellfun('isempty',strfind(RMS_scores.name, 'FORCEPLEIST'));

RMS_scores.DV_cosmo = nan*ones(length(ens_runs),1);
RMS_scores.MV_cosmo = nan*ones(length(ens_runs),1);
RMS_scores.LW_cosmo = nan*ones(length(ens_runs), 1);
RMS_scores.DH_cosmo = nan*ones(length(ens_runs), 1);

RMS_scores.DV_algae = nan*ones(length(ens_runs),1);
RMS_scores.MV_algae = nan*ones(length(ens_runs),1);
RMS_scores.LW_algae = nan*ones(length(ens_runs), 1);
RMS_scores.DH_algae = nan*ones(length(ens_runs), 1);

RMS_scores.total = nan*ones(length(ens_runs), 1);

%% Concatenate ages together for scoring
% For erratics Be-10 and bedrock C-14
DAN_UM_DV_BV_cosmo_ages = [data.parsed_output.DAN.erratics.t10St; ...
    data.parsed_output.UM.erratics.t10St; ...
    data.parsed_output.DV.erratics.t10St; ...
    data.parsed_output.BV.erratics.t10St];

DAN_UM_DV_BV_ages_index = DAN_UM_DV_BV_cosmo_ages<=-t_P2(1);

DAN_UM_DV_BV_cosmo_ages = DAN_UM_DV_BV_cosmo_ages(DAN_UM_DV_BV_ages_index);

DAN_UM_DV_BV_cosmo_elev = [DAN_cosmo2center; UM_cosmo2center; ...
    DV_cosmo2center; BV_cosmo2center];
DAN_UM_DV_BV_cosmo_elev = DAN_UM_DV_BV_cosmo_elev(DAN_UM_DV_BV_ages_index);

MV_cosmo_ages = [data.parsed_output.MV_walls.erratics.t10St; ...
    data.parsed_output.MV_floor.erratics.t10St];

MV_cosmo_ages_index = MV_cosmo_ages<=-t_P2(1);

MV_cosmo_elev = [MV_walls_cosmo2center; MV_floor_cosmo2center];
MV_cosmo_elev = MV_cosmo_elev(MV_cosmo_ages_index);
MV_cosmo_ages = MV_cosmo_ages(MV_cosmo_ages_index);

LW_cosmo_ages = [data.parsed_output.LW.erratics.t10St(LW_walls_cosmo_index); ...
    data.parsed_output.LW.erratics.t10St(~LW_walls_cosmo_index)]; 
LW_cosmo_ages_index = (LW_cosmo_ages<=-t_P2(1));

LW_cosmo_elev = [LW_walls_cosmo2center; LW_valley_cosmo2center];
LW_cosmo_ages = LW_cosmo_ages(LW_cosmo_ages_index);
LW_cosmo_elev = LW_cosmo_elev(LW_cosmo_ages_index);

DH_cosmo_ages = [data.parsed_output.DH.erratics.t10St; ...
    data.parsed_output.DH.bedrock.t14St];
DH_cosmo_ages_index = (DH_cosmo_ages<=-t_P2(1));

DH_cosmo_elev = [data.parsed_output.DH.erratics.HeightAboveIceMargin; ...
    data.parsed_output.DH.bedrock.HeightAboveIceMargin];
DH_cosmo_elev = DH_cosmo_elev(DH_cosmo_ages_index);
DH_cosmo_ages = DH_cosmo_ages(DH_cosmo_ages_index);

% For algae radiocarbon

DAN_UM_DV_BV_algae_ages = [DAN_algae.calYrBP; DVBV_algae.calYrBP];
DAN_UM_DV_BV_algae_elev = [DAN_algae2center; DVBV_algae2center];

MV_algae_ages = [MV_walls_algae.calYrBP; MV_floor_algae.calYrBP];
MV_algae_elev = [MV_walls_algae2center; MV_floor_algae2center];

LW_algae_ages = LW_algae.calYrBP;
LW_algae_elev = LW_algae2center;

DH_algae_ages = [DH_Darwin_algae.AveCalYr; DH_Diamond_algae.AveCalYr];
DH_algae_elev = [DH_Darwin_algae.HeightAboveIceMargin; ...
    DH_Diamond_algae.HeightAboveIceMargin];
    

%% Now loop through each run in the ensemble and score it against the geochronologic data,
%which includes the modern glacier surface
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
    plot(-t_P2, S_P2(:, DV_x_P_ind), plotline{:}, 'linewidth', 1.25);
    
    figure(2)
    hold on
    MV_ax = gca;
    plot(-t_P2, S_P2(:, MV_x_P_ind), plotline{:}, 'linewidth', 1.25);
    
    figure(3)
    hold on
    LW_ax = gca;
    plot(-t_P2, S_P2(:, LW_x_P_ind), plotline{:}, 'linewidth', 1.25);

    figure(4)
    hold on
    DH_ax = gca;
    plot(-t_P2, S_at_GL-S_at_GL(end), plotline{:}, 'linewidth', 1.25);
   
%% Now score this model run against geochronology

model_DV_cosmo = interp1(-t_P2, S_P2(:, DV_x_P_ind), ...
    DAN_UM_DV_BV_cosmo_ages);
model_DV_algae = interp1(-t_P2, S_P2(:, DV_x_P_ind), ...
    DAN_UM_DV_BV_algae_ages);

RMS_scores.DV_cosmo(jj) = sqrt(sum((DAN_UM_DV_BV_cosmo_elev - ... 
    model_DV_cosmo).^2)./length(DAN_UM_DV_BV_cosmo_elev));
RMS_scores.DV_algae(jj) = sqrt(sum((DAN_UM_DV_BV_algae_elev - ...
    model_DV_algae).^2)./length(DAN_UM_DV_BV_algae_elev));

model_MV_cosmo = interp1(-t_P2, S_P2(:, MV_x_P_ind), ... 
    MV_cosmo_ages);
model_MV_algae = interp1(-t_P2, S_P2(:, MV_x_P_ind), ... 
    MV_algae_ages);

RMS_scores.MV_cosmo(jj) = sqrt(sum((MV_cosmo_elev - ... 
    model_MV_cosmo).^2)./length(MV_cosmo_elev));
RMS_scores.MV_algae(jj) = sqrt(sum((MV_algae_elev - ... 
    model_MV_algae).^2)./length(MV_algae_elev));

model_LW_cosmo = interp1(-t_P2, S_P2(:, LW_x_P_ind), ...
    LW_cosmo_ages);
model_LW_algae = interp1(-t_P2, S_P2(:, LW_x_P_ind), ...
    LW_algae_ages);

RMS_scores.LW_cosmo(jj) = sqrt(sum((LW_cosmo_elev - ... 
    model_LW_cosmo).^2)./length(LW_cosmo_elev));
RMS_scores.LW_algae(jj) = sqrt(sum((LW_algae_elev - ... 
    model_LW_algae).^2)./length(LW_algae_elev));

model_DH_cosmo = interp1(-t_P, S_P(:, 1)-S_at_GL(end), ...
    DH_cosmo_ages);
model_DH_algae = interp1(-t_P, S_P(:, 1)-S_at_GL(end), ...
    DH_algae_ages);

RMS_scores.DH_cosmo(jj) = sqrt(sum((DH_cosmo_elev - ...
    model_DH_cosmo).^2)./length(DH_cosmo_elev));
RMS_scores.DH_algae(jj) = sqrt(sum((DH_algae_elev - ...
    model_DH_algae).^2)./length(DH_algae_elev));

% RMS_scores.total(jj) = (RMS_scores.DAN(jj).*RMS_scores.DV(jj).*RMS_scores.MV(jj)...
%     .*RMS_scores.LW(jj).*RMS_scores.DH(jj)).^(1/5);

if strfind(ens_runs{jj}, 'smooth9ka')
    RMS_scores.LR04_sealev(jj) = NaN;
end

end
%% total model score following section 2.4.1 (approach [a]) from Pollard et al., 2016

DV50_cosmo = median(RMS_scores.DV_cosmo(1:end-1));
MV50_cosmo = median(RMS_scores.MV_cosmo(1:end-1));
LW50_cosmo = median(RMS_scores.LW_cosmo(1:end-1));
DH50_cosmo = median(RMS_scores.DH_cosmo(1:end-1));

DV50_algae = median(RMS_scores.DV_algae(1:end-1));
MV50_algae = median(RMS_scores.MV_algae(1:end-1));
LW50_algae = median(RMS_scores.LW_algae(1:end-1));
DH50_algae = median(RMS_scores.DH_algae(1:end-1));

RMS_scores.cosmo_total = exp(-(RMS_scores.DV_cosmo./DV50_cosmo +...
    RMS_scores.MV_cosmo./MV50_cosmo + RMS_scores.LW_cosmo./LW50_cosmo ...
    + RMS_scores.DH_cosmo./DH50_cosmo));

RMS_scores.algae_total = exp(-(RMS_scores.DV_algae./DV50_algae +...
    RMS_scores.MV_algae./MV50_algae + RMS_scores.LW_algae./LW50_algae + ...
    RMS_scores.DH_algae./DH50_algae));
%% model score statistics for each parameter
tau1000_mean_cosmo = mean(RMS_scores.cosmo_total(RMS_scores.tau==1000));
tau1000_mean_algae = mean(RMS_scores.algae_total(RMS_scores.tau==1000));
tau1000_std_cosmo = std(RMS_scores.cosmo_total(RMS_scores.tau==1000));
tau1000_std_algae = std(RMS_scores.algae_total(RMS_scores.tau==1000));
tau2000_mean_cosmo = mean(RMS_scores.cosmo_total(RMS_scores.tau==2000));
tau2000_mean_algae = mean(RMS_scores.algae_total(RMS_scores.tau==2000));
tau2000_std_cosmo = std(RMS_scores.cosmo_total(RMS_scores.tau==2000));
tau2000_std_algae = std(RMS_scores.algae_total(RMS_scores.tau==2000));

ocfac0_5_mean_cosmo = mean(RMS_scores.cosmo_total(RMS_scores.ocfac==0.5));
ocfac0_5_std_cosmo = std(RMS_scores.cosmo_total(RMS_scores.ocfac==0.5));
ocfac1_mean_cosmo = mean(RMS_scores.cosmo_total(RMS_scores.ocfac==1));
ocfac1_std_cosmo = std(RMS_scores.cosmo_total(RMS_scores.ocfac==1));
ocfac0_5_mean_algae = mean(RMS_scores.algae_total(RMS_scores.ocfac==0.5));
ocfac0_5_std_algae = std(RMS_scores.algae_total(RMS_scores.ocfac==0.5));
ocfac1_mean_algae = mean(RMS_scores.algae_total(RMS_scores.ocfac==1));
ocfac1_std_algae = std(RMS_scores.algae_total(RMS_scores.ocfac==1));


shelf5_mean_cosmo = mean(RMS_scores.cosmo_total(RMS_scores.crhshelf == 1e-5));
shelf5_std_cosmo = std(RMS_scores.cosmo_total(RMS_scores.crhshelf == 1e-5));
shelf6_mean_cosmo = mean(RMS_scores.cosmo_total(RMS_scores.crhshelf == 1e-6));
shelf6_std_cosmo = std(RMS_scores.cosmo_total(RMS_scores.crhshelf == 1e-6));
shelf7_mean_cosmo = mean(RMS_scores.cosmo_total(RMS_scores.crhshelf == 1e-7));
shelf7_std_cosmo = std(RMS_scores.cosmo_total(RMS_scores.crhshelf == 1e-7));
shelf5_mean_algae = mean(RMS_scores.algae_total(RMS_scores.crhshelf == 1e-5));
shelf5_std_algae = std(RMS_scores.algae_total(RMS_scores.crhshelf == 1e-5));
shelf6_mean_algae = mean(RMS_scores.algae_total(RMS_scores.crhshelf == 1e-6));
shelf6_std_algae = std(RMS_scores.algae_total(RMS_scores.crhshelf == 1e-6));
shelf7_mean_algae = mean(RMS_scores.algae_total(RMS_scores.crhshelf == 1e-7));
shelf7_std_algae = std(RMS_scores.algae_total(RMS_scores.crhshelf == 1e-7));

calv0_7_mean_cosmo = mean(RMS_scores.cosmo_total(RMS_scores.calv == 0.7));
calv0_7_std_cosmo = std(RMS_scores.cosmo_total(RMS_scores.calv == 0.7));
calv1_mean_cosmo = mean(RMS_scores.cosmo_total(RMS_scores.calv == 1));
calv1_std_cosmo = std(RMS_scores.cosmo_total(RMS_scores.calv == 1));
calv0_7_mean_algae = mean(RMS_scores.algae_total(RMS_scores.calv == 0.7));
calv0_7_std_algae = std(RMS_scores.algae_total(RMS_scores.calv == 0.7));
calv1_mean_algae = mean(RMS_scores.algae_total(RMS_scores.calv == 1));
calv1_std_algae = std(RMS_scores.algae_total(RMS_scores.calv == 1));

FORCEPLEIST_mean_cosmo = mean(RMS_scores.cosmo_total(RMS_scores.LR04_sealev == true));
FORCEPLEIST_std_cosmo = std(RMS_scores.cosmo_total(RMS_scores.LR04_sealev == true));
SPRATT_mean_cosmo = mean(RMS_scores.cosmo_total(RMS_scores.LR04_sealev == false));
SPRATT_std_cosmo = std(RMS_scores.cosmo_total(RMS_scores.LR04_sealev == false));

FORCEPLEIST_mean_algae = mean(RMS_scores.algae_total(RMS_scores.LR04_sealev == true));
FORCEPLEIST_std_algae = std(RMS_scores.algae_total(RMS_scores.LR04_sealev == true));
SPRATT_mean_algae = mean(RMS_scores.algae_total(RMS_scores.LR04_sealev == false));
SPRATT_std_algae = std(RMS_scores.algae_total(RMS_scores.LR04_sealev == false));
%% plot extracted datapoints

figure(1)
    hold on
    title({'\it Dubris and Bibra Valleys,'; 'Danum Platform, Updog Mtn '}, 'FontWeight', 'normal')
    DAN_UM_DV_BV_cosmo_plot = plot(DAN_UM_DV_BV_cosmo_ages, DAN_UM_DV_BV_cosmo_elev, 'ko', 'MarkerFaceColor', [0.8, 0.8, 0.8]);
    DAN_UM_DV_BV_algae_plot = plot(DAN_UM_DV_BV_algae_ages, DAN_UM_DV_BV_algae_elev, 'ko', 'MarkerFaceColor', [0.2, 0.8, 0.8]);
    set(gca, 'Fontsize', 15)
    set(gcf, 'Position', [202   438   292   260], 'PaperPositionMode', 'auto')
    ylabel('Elevation (m)')
    xlabel ('Age (kyr BP)')
    legend([DAN_UM_DV_BV_cosmo_plot, DAN_UM_DV_BV_algae_plot], ...
        {'^1^0Be Exposure ages'; 'Algae ^1^4C ages'})
    grid on
    
figure(2)
    hold on
    title('\it Magnis Valley ', 'FontWeight', 'normal')
    plot(MV_cosmo_ages, MV_cosmo_elev, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8]);
    plot(MV_algae_ages, MV_algae_elev, 'ko', 'MarkerFaceColor', [0.2 0.8 0.8]);
    set(gca, 'Fontsize', 15)
    set(gcf, 'Position', [202   438   292   260], 'PaperPositionMode', 'auto')
    ylabel('Elevation (m)')
    xlabel ('Age (kyr BP)')
    grid on
    
figure(3)
    hold on
    title('\it Lake Wellman ', 'FontWeight', 'normal')
    plot(LW_cosmo_ages, LW_cosmo_elev, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8]);
    plot(LW_algae_ages, LW_algae_elev, 'ko', 'MarkerFaceColor', [0.2 0.8 0.8]);
    set(gca, 'Fontsize')
    set(gcf, 'Position', [202   438   292   260], 'PaperPositionMode', 'auto')
    ylabel('Elevation (m)')
    xlabel ('Age (kyr BP)')
    grid on
    
figure(4)
    hold on
    title('\it Diamond Hill ', 'FontWeight', 'normal')
    plot(DH_cosmo_ages, DH_cosmo_elev, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8])
    plot(DH_algae_ages, DH_algae_elev, 'ko', 'MarkerFaceColor', [0.2 0.8 0.8])
    set(gca, 'Fontsize')
    set(gcf, 'Position', [202   438   292   260], 'PaperPositionMode', 'auto')
    ylabel('Height above modern glacier (m)')
    xlabel ('Age (kyr BP)')
    grid on
    
%% Algae: get the five highest scores and plot them in a better color
highest_scores_algae = sort(RMS_scores.algae_total);
highest_scores_algae = highest_scores_algae(end-5:end);
highest_scores_cosmo = sort(RMS_scores.cosmo_total);
highest_scores_cosmo = highest_scores_cosmo(end-5:end);

% this is a nice purple and green colormap, but don't use it here.
% colors = [118,42,131
% 175,141,195
% 231,212,232
% 217,240,211
% 127,191,123
% 27,120,55]./255;

cosmo_markerstyle = {'o', 'd', 's', '<', '^', '>'};
mark_pos = 1; % marker position so they will be offset to avoid overlap.
algae_cmap = summer(length(highest_scores_algae));

for jj = 1:length(highest_scores_algae);
highest_scores_algae_index(jj) = find(RMS_scores.algae_total == highest_scores_algae(jj));
highest_scores_cosmo_index(jj) = find(RMS_scores.cosmo_total == highest_scores_cosmo(jj));

load([runs_dir,'/', ens_runs{highest_scores_algae_index(jj)}], 't_P', 't_P2', ...
        'x_P', 'x_P2', 'S_P', 'S_P2', 'S_at_GL');

      

    
    figure(1)

    hold on
    DV_ax = gca;
    plot(-t_P2, S_P2(:, DV_x_P_ind), 'Color', algae_cmap(jj,:), 'linewidth', 2);
    
    figure(2)
    hold on
    plot(-t_P2, S_P2(:, MV_x_P_ind), 'Color', algae_cmap(jj,:),  'linewidth', 2);
    
    figure(3)
    hold on
    LW_ax = gca;
    plot(-t_P2, S_P2(:, LW_x_P_ind), 'Color', algae_cmap(jj,:), 'linewidth', 2);

    figure(4)
    hold on
    DH_ax = gca;
    plot(-t_P2, S_at_GL-S_at_GL(end), 'Color', algae_cmap(jj,:),  'linewidth', 2);
 
load([runs_dir,'/', ens_runs{highest_scores_cosmo_index(jj)}], 't_P', 't_P2', ...
        'x_P', 'x_P2', 'S_P', 'S_P2', 'S_at_GL');
    
    figure(1)
    hold on
    DV_ax = gca;
    plot(-t_P2(mark_pos:10:end), S_P2(mark_pos:10:end, DV_x_P_ind), ...
        cosmo_markerstyle{jj}, 'color', 'k');
    
    figure(2)
    hold on
    plot(-t_P2(mark_pos:10:end), S_P2((mark_pos:10:end), MV_x_P_ind), ...
        cosmo_markerstyle{jj}, 'color', 'k');
    
    figure(3)
    hold on
    LW_ax = gca;
    plot(-t_P2(mark_pos:10:end), S_P2((mark_pos:10:end), LW_x_P_ind), ... 
        cosmo_markerstyle{jj}, 'color', 'k');

    figure(4)
    hold on
    DH_ax = gca;
    plot(-t_P2(mark_pos:10:end), S_at_GL(mark_pos:10:end)-S_at_GL(end), ...
        cosmo_markerstyle{jj}, 'color', 'k');
    
    mark_pos = mark_pos + 7; % offset markers for clarity
end

%% plot model score stats if the runs are from a PSU ice model ensemble
if ~isempty(strfind(runs_dir, 'sealev_ens')) 
    
    
figure(6)
subplot(1,4,1)
FaceColor = [188,189,220]./255;
bar([1, 0], [mean(RMS_scores.algae_total(RMS_scores.LR04_sealev == 1)),...
    mean(RMS_scores.algae_total(RMS_scores.LR04_sealev == 0))], 'FaceColor', FaceColor)
set(gca, 'XTickLabel', {'Spratt', 'LR04'}, 'Fontsize', 15)
xlabel('Sea-level forcing')
ylim([0 0.035])
xlim([-0.5 1.5])
ylabel('Mean model score')
grid on

subplot(1,4,2)
bar([1000, 2000], [mean(RMS_scores.algae_total(RMS_scores.tau==1000)),...
    mean(RMS_scores.algae_total(RMS_scores.tau==2000))], 'FaceColor', FaceColor)
set(gca, 'fontsize', 15, 'YTickLabel', {'' '' '' ''})
xlabel({'\tau: rebound'; 'response time (yr)'})
ylim([0 0.035])
xlim([500 2500])
grid on

subplot(1,4,3)
bar([0.5 1], [mean(RMS_scores.algae_total(RMS_scores.ocfac==0.5))...
    mean(RMS_scores.algae_total(RMS_scores.ocfac==1))], 'FaceColor', FaceColor); 
xlabel({'Sub-shelf'; 'melting factor'})
set(gca, 'fontsize', 15, 'Xtick', [0.5 1], 'YTickLabel', {'' '' '' ''})
ylim([0 0.035])
xlim([0.25 1.25])
grid on

subplot(1,4,4)
 bar([-5 -6 -7], [mean(RMS_scores.algae_total(RMS_scores.crhshelf == 1e-5))...
     mean(RMS_scores.algae_total(RMS_scores.crhshelf == 1e-6))...
     mean(RMS_scores.algae_total(RMS_scores.crhshelf == 1e-7))], 'FaceColor', FaceColor)
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

highest_scoring_runs = regexprep(RMS_scores.name(highest_scores_algae_index),...
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
end