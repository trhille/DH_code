% This script prepares geochronology data for evaluating output of
% flow-band model runs. 
load cosmo_data.mat
load algae_data.mat


%% DV, BV, DAN, UM cosmo ages. Only include retreat ages in the polynomial fit. Limit ages will be fit
% separately
n_DVBV_cosmo = 1;
DVBV_cosmo_retreat_start = 7.8e3;
DVBV_cosmo_poly = polyfit([data.parsed_output.UM.erratics.t10St(data.parsed_output.UM.erratics.t10St<=DVBV_cosmo_retreat_start);...
    data.parsed_output.DV.erratics.t10St(data.parsed_output.DV.erratics.t10St<=DVBV_cosmo_retreat_start);...
    data.parsed_output.BV.erratics.t10St(data.parsed_output.BV.erratics.t10St<=DVBV_cosmo_retreat_start);...
    data.parsed_output.DAN.erratics.t10St(data.parsed_output.DAN.erratics.t10St<=DVBV_cosmo_retreat_start)], ...
    [data.parsed_output.UM.erratics.elev(data.parsed_output.UM.erratics.t10St<=DVBV_cosmo_retreat_start);...
    data.parsed_output.DV.erratics.elev(data.parsed_output.DV.erratics.t10St<=DVBV_cosmo_retreat_start); ...
    data.parsed_output.BV.erratics.elev(data.parsed_output.BV.erratics.t10St<=DVBV_cosmo_retreat_start);...
    data.parsed_output.DAN.erratics.elev(data.parsed_output.DAN.erratics.t10St<=DVBV_cosmo_retreat_start)], n_DVBV_cosmo);

DVBV_cosmo_time = 3e3:100:DVBV_cosmo_retreat_start;
DVBV_cosmo_fit = polyval(DVBV_cosmo_poly, DVBV_cosmo_time);
DVBV_cosmo_rate = gradient(DVBV_cosmo_fit)./gradient(DVBV_cosmo_time);

figure(1); clf
hold on; grid on; box on
herrorbar(data.parsed_output.BV.erratics.t10St,data.parsed_output.BV.erratics.elev,...
    data.parsed_output.BV.erratics.int10St, 'bo')
herrorbar(data.parsed_output.DV.erratics.t10St,data.parsed_output.DV.erratics.elev,...
    data.parsed_output.DV.erratics.int10St, 'bo')
herrorbar(data.parsed_output.UM.erratics.t10St,data.parsed_output.UM.erratics.elev,...
    data.parsed_output.UM.erratics.int10St, 'ko')
herrorbar(data.parsed_output.DAN.erratics.t10St,data.parsed_output.DAN.erratics.elev,...
    data.parsed_output.DAN.erratics.int10St, 'ko')

plot(DVBV_cosmo_time, DVBV_cosmo_fit);
xlim([0 1e4])
xlabel ('Exposure Age (yr BP)')
ylabel ('Elevation (m)')
set(gca, 'FontSize', 14)

DAN_UM_linear = slmengine([data.parsed_output.DAN.erratics.t10St;...
    data.parsed_output.UM.erratics.t10St],...
    [data.parsed_output.DAN.erratics.elev;...
    data.parsed_output.UM.erratics.elev],'degree',...
    'linear','knots',[5000, 7000, 15000], 'plot','off');

DAN_UM_linear_fit_time = sort([data.parsed_output.DAN.erratics.t10St;...
    data.parsed_output.UM.erratics.t10St]);
DAN_UM_linear_fit = slmeval(DAN_UM_linear_fit_time, DAN_UM_linear);

DVBV_cosmo_linear = slmengine([data.parsed_output.DV.erratics.t10St(data.parsed_output.DV.erratics.t10St<=1e4);...
    data.parsed_output.BV.erratics.t10St(data.parsed_output.BV.erratics.t10St<=1e4)],...
    [data.parsed_output.DV.erratics.elev(data.parsed_output.DV.erratics.t10St<=1e4);...
    data.parsed_output.BV.erratics.elev(data.parsed_output.BV.erratics.t10St<=1e4)],'degree',...
    'linear','knots',[2000, 7500, 10000], 'plot','off');

DVBV_cosmo_linear_fit_time = sort...
    ([data.parsed_output.DV.erratics.t10St(data.parsed_output.DV.erratics.t10St<=1e4);...
    data.parsed_output.BV.erratics.t10St(data.parsed_output.BV.erratics.t10St<=1e4)]);

DVBV_cosmo_linear_fit = slmeval(DVBV_cosmo_linear_fit_time,...
    DVBV_cosmo_linear);

plot(DVBV_cosmo_linear_fit_time,...
    DVBV_cosmo_linear_fit,'k', 'linewidth', 2);
plot(DAN_UM_linear_fit_time, DAN_UM_linear_fit, 'r', 'linewidth', 2);
% S_at_GL_linear = [slmeval(time, slm_linear), max(all_height)]; 
% 
%% DV, DAN, and BV algae alges
n_DVBV_algae = 1;
DVBV_algae_retreat_start = 1e4;
DVBV_algae_time = 0:100:DVBV_algae_retreat_start;
DVBV_algae_poly = polyfit([DVBV_algae.calYrBP(DVBV_algae.calYrBP<=DVBV_algae_retreat_start);...
    DAN_algae.calYrBP(DAN_algae.calYrBP<=DVBV_algae_retreat_start)],...
    [DVBV_algae.Elevation(DVBV_algae.calYrBP<=DVBV_algae_retreat_start);...
    DAN_algae.Elevation(DAN_algae.calYrBP<=DVBV_algae_retreat_start)], n_DVBV_algae);

DVBV_algae_fit =  polyval(DVBV_algae_poly,DVBV_algae_time);
DVBV_algae_rate = gradient(DVBV_algae_fit)./gradient(DVBV_algae_time);

figure(2); clf
hold on; grid on; box on
plot(DAN_algae.calYrBP, DAN_algae.Elevation, 'ks')
plot(DVBV_algae.calYrBP, DVBV_algae.Elevation, 'bs')
plot(DVBV_algae_time, DVBV_algae_fit)
xlabel ('Age (cal yr BP)')
ylabel ('Elevation (m)')
set(gca, 'FontSize', 14)

DVBV_algae_linear = slmengine([DAN_algae.calYrBP; DVBV_algae.calYrBP],...
    [DAN_algae.Elevation; DVBV_algae.Elevation],'degree',...
    'linear','knots',[0, 5500, 9000, 15000], 'plot','off');

DVBV_algae_linear_fit_time = sort([DAN_algae.calYrBP; DVBV_algae.calYrBP]);
DVBV_algae_linear_fit = slmeval(DVBV_algae_linear_fit_time,...
    DVBV_algae_linear);

plot(DVBV_algae_linear_fit_time, DVBV_algae_linear_fit,'k', 'linewidth', 2);

%% compare algae and cosmo rates in Dubris-Bibra:
figure(3)
clf 
hold on; grid on; box on
plot(DVBV_algae_time, DVBV_algae_rate, 'k--', 'linewidth', 2)
plot(DVBV_cosmo_time, DVBV_cosmo_rate, 'r--', 'linewidth', 2)
legend ('Algae thinning rate', 'Cosmo thinning rate')
xlabel ('Time (yr BP)')
ylabel ('Thinning rate (m/yr)')
set(gca, 'FontSize', 14)


%% Magnis Valley exposure ages
% MV_floor and MV_walls must be fit separately
n_MV_floor_cosmo = 1;
MV_floor_retreat_start = 1e4; 
figure(4)
clf; hold on; box on; grid on;
herrorbar(data.parsed_output.MV_floor.erratics.t10St,...
    data.parsed_output.MV_floor.erratics.elev,...
    data.parsed_output.MV_floor.erratics.int10St, 'bo')
herrorbar(data.parsed_output.MV_walls.erratics.t10St,...
     data.parsed_output.MV_walls.erratics.elev,...
    data.parsed_output.MV_walls.erratics.int10St, 'ko')


MV_floor_cosmo_poly = polyfit(data.parsed_output.MV_floor.erratics.t10St...
    (data.parsed_output.MV_floor.erratics.t10St<=MV_floor_retreat_start),...
    data.parsed_output.MV_floor.erratics.elev...
    (data.parsed_output.MV_floor.erratics.t10St<=MV_floor_retreat_start),...
    n_MV_floor_cosmo);

MV_floor_cosmo_time = 5e3:100:1e4;
MV_floor_cosmo_fit = polyval(MV_floor_cosmo_poly, MV_floor_cosmo_time);



%MV_walls
n_MV_walls_cosmo = 1;
MV_walls_retreat_start = 1e4; 


MV_walls_cosmo_poly = polyfit(data.parsed_output.MV_walls.erratics.t10St...
    (data.parsed_output.MV_walls.erratics.t10St<=MV_walls_retreat_start),...
    data.parsed_output.MV_walls.erratics.elev...
    (data.parsed_output.MV_walls.erratics.t10St<=MV_walls_retreat_start),...
    n_MV_walls_cosmo);

MV_walls_cosmo_time = 5e3:100:MV_walls_retreat_start;
MV_walls_cosmo_fit = polyval(MV_walls_cosmo_poly, MV_walls_cosmo_time);

figure(4)
plot(MV_floor_cosmo_time, MV_floor_cosmo_fit)
plot(MV_walls_cosmo_time, MV_walls_cosmo_fit)

MV_floor_outlier = {'14-HAT-059-MV', '14-HAT-045-MV'};
MV_floor_outlier_index = ismember(data.parsed_output.MV_floor.erratics.SampleName,...
    MV_floor_outlier);

MV_floor_cosmo_linear = slmengine(...
    data.parsed_output.MV_floor.erratics.t10St(~MV_floor_outlier_index),...
    data.parsed_output.MV_floor.erratics.elev(~MV_floor_outlier_index),'degree',...
    'linear','knots',[5e3, 9500, 15000], 'plot','off');

MV_walls_cosmo_linear = slmengine(...
    data.parsed_output.MV_walls.erratics.t10St,...
    data.parsed_output.MV_walls.erratics.elev,'degree',...
    'linear','knots',[5e3, 8e3, 15000], 'plot','off');

MV_floor_cosmo_linear_fit_time = ...
    sort(data.parsed_output.MV_floor.erratics.t10St(~MV_floor_outlier_index));
MV_walls_cosmo_linear_fit_time = ...
    sort(data.parsed_output.MV_walls.erratics.t10St);

MV_floor_cosmo_linear_fit = slmeval(MV_floor_cosmo_linear_fit_time,...
    MV_floor_cosmo_linear);
MV_walls_cosmo_linear_fit = slmeval(MV_walls_cosmo_linear_fit_time,...
    MV_walls_cosmo_linear);

plot(MV_floor_cosmo_linear_fit_time, MV_floor_cosmo_linear_fit)
plot(MV_walls_cosmo_linear_fit_time, MV_walls_cosmo_linear_fit)

%% MV algae ages
figure(4)
plot(MV_floor_algae.calYrBP, MV_floor_algae.Elevation, 'rs')
MV_floor_algae_linear = slmengine(...
    MV_floor_algae.calYrBP(MV_floor_algae.calYrBP<=1e4),...
    MV_floor_algae.Elevation(MV_floor_algae.calYrBP<=1e4),...
    'degree',...
    'linear','knots',[1e3,6750, 8e3, 15000], 'plot','off');
MV_floor_algae_linear_fit_time = ...
    sort(MV_floor_algae.calYrBP(MV_floor_algae.calYrBP<=1e4));

MV_floor_algae_linear_fit = slmeval(MV_floor_algae_linear_fit_time,...
    MV_floor_algae_linear);

plot(MV_floor_algae_linear_fit_time, MV_floor_algae_linear_fit)


%% LW exposure ages
%call everything >12.5 kyr pre-exposed

LW_cosmo_use_ind = logical(1-(data.parsed_output.LW.erratics.t10St>=12.5e3));

LW_cosmo_linear = slmengine(...
    data.parsed_output.LW.erratics.t10St(LW_cosmo_use_ind),...
    data.parsed_output.LW.erratics.elev(LW_cosmo_use_ind),'degree',...
    'linear','knots',[0e3, 2e3, 8e3, 15000], 'plot','off');

LW_cosmo_linear_fit_time = ...
    sort(data.parsed_output.LW.erratics.t10St(LW_cosmo_use_ind));
LW_cosmo_linear_fit = slmeval(LW_cosmo_linear_fit_time, LW_cosmo_linear);

figure(5); clf
hold on; box on; grid on
herrorbar(data.parsed_output.LW.erratics.t10St, ...
    data.parsed_output.LW.erratics.elev,...
    data.parsed_output.LW.erratics.int10St, 'ko')
plot(LW_cosmo_linear_fit_time, LW_cosmo_linear_fit)

%% LW algae radiocarbon

figure(5); hold on
plot(LW_algae.calYrBP, LW_algae.Elevation, 'rs')

LW_outliers =...
{'LW-13-164A';
'LW-13-164B';
'LW-13-168';
'LW-13-170'}; %samples from Darwin Depression whose elevations are too low

LW_algae_use_ind = ~ismember(LW_algae.SampleID, LW_outliers);
LW_algae_linear = slmengine(...
    LW_algae.calYrBP(LW_algae_use_ind),...
    LW_algae.Elevation(LW_algae_use_ind),'degree',...
    'linear','knots',[2e3, 8.5e3, 15000], 'plot','off');

LW_algae_linear_fit_time = sort(LW_algae.calYrBP);
LW_algae_linear_fit = slmeval(LW_algae_linear_fit_time, LW_algae_linear);

plot(LW_algae_linear_fit_time, LW_algae_linear_fit, '-')

%% LGM elevation estimate from a linear projection between limits in valleys and on walls

close all
project_samples_to_centerline

%% Save constraints
% For each fit, subtract the maximum elevation to get thinning relative to
% LGM. This will be the actual metric you should use to score the runs
save geochron_constraints.mat *linear_fit *linear_fit_time -append









