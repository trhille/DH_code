% This script projects algae and exposure ages to the glacier centerline by
% fitting a parabolic cross-glacier surface profile to LGM limits in
% valleys and on valley walls. 

load algae_data.mat
load cosmo_data.mat

addpath('/Users/trevorhillebrand/Documents/MATLAB/Toolboxes/AntarcticMappingTools_v5.00//AntarcticMappingTools');

centerline = shaperead(['/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton/Modeling/',...
    'Koutnik model/DH_code/DH_DATA/Shapefiles/Hatherton Flowline/Hatherton flowline.shp']);
%interpolate the centerline, since it is only defined at 18 points in the
%shapefile
centerline_x_interp = interp1((1:length(centerline.X)),...
    centerline.X, linspace(1, length(centerline.X), 1000));
centerline_y_interp = interp1((1:length(centerline.Y)),...
    centerline.Y, linspace(1, length(centerline.Y), 1000));



%% Lake Wellman

%define which samples consitute the LGM limit in the valley and on the
%walls nearest the glacier
LW_valley_limit_cosmo = {'13-HAT-030-LW', '13-HAT-029-LW', '13-HAT-031-LW'};
LW_walls_limit_cosmo = {'13-HAT-044-LW', '13-HAT-047-LW', '13-HAT-006-LW',...
    '13-HAT-007-LW'};

LW_walls_limit_cosmo_index = ismember(data.parsed_output.LW.erratics.SampleName,...
    LW_walls_limit_cosmo);
LW_valley_limit_cosmo_index = ismember(data.parsed_output.LW.erratics.SampleName,...
    LW_valley_limit_cosmo);

%indentify all samples that count as being on valley walls. This is
%somewhat ambiguous.
LW_walls_cosmo_samples = {'13-HAT-004-LW'; '13-HAT-006-LW'; '13-HAT-007-LW';...
    '13-HAT-036-LW'; '13-HAT-041-LW'; '13-HAT-042-LW'; '13-HAT-044-LW'; 
    '13-HAT-047-LW'};
LW_walls_cosmo_index = ismember(data.parsed_output.LW.erratics.SampleName,...
    LW_walls_cosmo_samples);

%Index into the big data structure to access only those samples
LW_valley_index = ismember(data.output.SampleName,...
    LW_valley_limit_cosmo);
LW_walls_index = ismember(data.output.SampleName,...
    LW_walls_limit_cosmo);

%find the elevations...
LW_valley_limit_elev = data.output.elev(LW_valley_index);
LW_walls_limit_elev = data.output.elev(LW_walls_index);

%...and lat-lons
LW_valley_limit_lat = data.input.lat(ismember(data.input.SampleName,...
    LW_valley_limit_cosmo));
LW_walls_limit_lat = data.input.lat(ismember(data.input.SampleName,...
    LW_walls_limit_cosmo));

LW_valley_limit_long = data.input.long(ismember(data.input.SampleName,...
    LW_valley_limit_cosmo));
LW_walls_limit_long = data.input.long(ismember(data.input.SampleName,...
    LW_walls_limit_cosmo));

% Convert to polar stereographic for calculation measurement
[LW_valley_limit_x, LW_valley_limit_y] = ll2ps(LW_valley_limit_lat, ...
    LW_valley_limit_long);

[LW_walls_limit_x, LW_walls_limit_y] = ll2ps(LW_walls_limit_lat, ...
    LW_walls_limit_long);


%calculate euclidean distance to the Hatherton glacier centerline for each
%sample. This is probably a good enough approximation
for jj = 1:length(LW_valley_limit_cosmo)
    LW_valley_distance_to_centerline(jj) = min(sqrt((centerline_x_interp - ...
        LW_valley_limit_x(jj)).^2 + (centerline_y_interp - LW_valley_limit_y(jj)).^2));
    
end


for kk = 1:length(LW_walls_limit_cosmo)
    LW_walls_distance_to_centerline(kk) = min(sqrt((centerline_x_interp - ...
        LW_walls_limit_x(kk)).^2 + (centerline_y_interp - LW_walls_limit_y(kk)).^2));
   
end

%Concatenate all together, since we are fitting one curve to these two
%populations

LW_limit_distance_to_centerline = [LW_valley_distance_to_centerline, LW_walls_distance_to_centerline];
LW_limit_elev = [LW_valley_limit_elev', LW_walls_limit_elev'];

%plot it up to make sure it looks reasonable
figure(1)
plot(LW_limit_distance_to_centerline, LW_limit_elev, 'ko'); hold on

% Good! Now fit a parabola to this, and evaluate it out to estimate the
% centerline thickness. We mirror elevations across the glacier to ensure
% a concave-down parabola, symmetric about the glacier centerline. Probably
% a little heavy-handed, but we are already making a biggish assumption
% that the cross surface profile is a parabola and that the glacier width
% does not change through time (ie., centerline stays put!).
n = 1; %2nd order polynomial
if n == 2
    LW_poly = polyfit([-LW_limit_distance_to_centerline,...
        LW_limit_distance_to_centerline] , [LW_limit_elev, LW_limit_elev],n);
elseif n == 1
    LW_poly = polyfit(LW_limit_distance_to_centerline, LW_limit_elev, n);

end
 LW_surface_profile = polyval(LW_poly, linspace(0, max(LW_limit_distance_to_centerline), 100));

LGM_surface_LW = LW_surface_profile(1);
% now plot up the surface profile to make sure it's a good fit
plot(linspace(0, max(LW_limit_distance_to_centerline), 100), LW_surface_profile)
title ('Lake Wellman')
xlabel ('Distance to centerline (m)')
ylabel ('Elevation (m)')
set (gca, 'fontsize', 13)

%% Magnis Valley 
%(copy-pasted from Lake Wellman and then edited)

%define which samples consitute the LGM limit in the valley and on the
%walls nearest the glacier
MV_valley_limit_cosmo = {'14-HAT-046-MV', '14-HAT-045-MV', '14-HAT-040-MV',...
    '14-HAT-041-MV', '14-HAT-042-MV', '14-HAT-044-MV'};
MV_walls_limit_cosmo = {'14-HAT-070-MV', '14-HAT-048-MV', '14-HAT-060-MV'};
   
%Index into the big data structure to access only those samples
MV_valley_index = ismember(data.output.SampleName,...
    MV_valley_limit_cosmo);
MV_walls_index = ismember(data.output.SampleName,...
    MV_walls_limit_cosmo);

%find the elevations...
MV_valley_limit_elev = data.output.elev(MV_valley_index);
MV_walls_limit_elev = data.output.elev(MV_walls_index);

%...and lat-lons
MV_valley_limit_lat = data.input.lat(ismember(data.input.SampleName,...
    MV_valley_limit_cosmo));
MV_walls_limit_lat = data.input.lat(ismember(data.input.SampleName,...
    MV_walls_limit_cosmo));

MV_valley_limit_long = data.input.long(ismember(data.input.SampleName,...
    MV_valley_limit_cosmo));
MV_walls_limit_long = data.input.long(ismember(data.input.SampleName,...
    MV_walls_limit_cosmo));

% Convert to polar stereographic for calculation measurement
[MV_valley_limit_x, MV_valley_limit_y] = ll2ps(MV_valley_limit_lat, ...
    MV_valley_limit_long);

[MV_walls_limit_x, MV_walls_limit_y] = ll2ps(MV_walls_limit_lat, ...
    MV_walls_limit_long);


%calculate euclidean distance to the Hatherton glacier centerline for each
%sample. This is probably a good enough approximation
for jj = 1:length(MV_valley_limit_cosmo)
    MV_valley_distance_to_centerline(jj) = min(sqrt((centerline_x_interp - ...
        MV_valley_limit_x(jj)).^2 + (centerline_y_interp - MV_valley_limit_y(jj)).^2));
    
end


for kk = 1:length(MV_walls_limit_cosmo)
    MV_walls_distance_to_centerline(kk) = min(sqrt((centerline_x_interp - ...
        MV_walls_limit_x(kk)).^2 + (centerline_y_interp - MV_walls_limit_y(kk)).^2));
   
end

%Concatenate all together, since we are fitting one curve to these two
%populations

MV_limit_distance_to_centerline = [MV_valley_distance_to_centerline,...
    MV_walls_distance_to_centerline];
MV_limit_elev = [MV_valley_limit_elev', MV_walls_limit_elev'];

%plot it up to make sure it looks reasonable
figure(2)
plot(MV_limit_distance_to_centerline, MV_limit_elev, 'ko'); hold on

% Good! Now fit a parabola to this, and evaluate it out to estimate the
% centerline thickness. We mirror elevations across the glacier to ensure
% a concave-down parabola, symmetric about the glacier centerline. Probably
% a little heavy-handed, but we are already making a biggish assumption
% that the cross surface profile is a parabola and that the glacier width
% does not change through time (ie., centerline stays put!).
n = 1; %2nd order polynomial
if n == 2
    MV_poly = polyfit([-MV_limit_distance_to_centerline,...
        MV_limit_distance_to_centerline] , [MV_limit_elev, MV_limit_elev],n);
elseif n == 1
	MV_poly = polyfit(MV_limit_distance_to_centerline,...
        MV_limit_elev,n); 
end

MV_surface_profile = polyval(MV_poly, linspace(0,...
    max(MV_limit_distance_to_centerline), 100));

LGM_surface_MV = MV_surface_profile(1);
% now plot up the surface profile to make sure it's a good fit
plot(linspace(0, max(MV_limit_distance_to_centerline), 100), MV_surface_profile)
title ('Magnis Valley')
xlabel ('Distance to centerline (m)')
ylabel ('Elevation (m)')
set (gca, 'fontsize', 13)

%% Dubris Valley 
%(copy-pasted from Lake Wellman and then edited)

%define which samples consitute the LGM limit in the valley and on the
%walls nearest the glacier
DV_limit_cosmo = {'13-HAT-067-DV', '13-HAT-068-DV', '13-HAT-069-DV'};
UM_limit_cosmo = {'13-HAT-137-UM', '13-HAT-138-UM', '13-HAT-140-UM'};

   
%Index into the big data structure to access only those samples
DV_limit_index = ismember(data.output.SampleName,...
    DV_limit_cosmo);
UM_limit_index = ismember(data.output.SampleName,...
    UM_limit_cosmo);

%the only way to find a substring in cell arrays. Very dumb:
DV_index = find(not(cellfun('isempty',strfind(data.output.SampleName, 'DV'))));
BV_index = find(not(cellfun('isempty',strfind(data.output.SampleName, 'BV'))));
DAN_index = find(not(cellfun('isempty',strfind(data.output.SampleName, 'DAN'))));
UM_index = find(not(cellfun('isempty',strfind(data.output.SampleName, 'UM'))));

%find the elevations...
DV_limit_elev = data.output.elev(DV_limit_index);
UM_limit_elev = data.output.elev(UM_limit_index);

%...and lat-lons
DV_limit_lat = data.input.lat(ismember(data.input.SampleName,...
    DV_limit_cosmo));
UM_limit_lat = data.input.lat(ismember(data.input.SampleName,...
    UM_limit_cosmo));

DV_limit_long = data.input.long(ismember(data.input.SampleName,...
    DV_limit_cosmo));
UM_limit_long = data.input.long(ismember(data.input.SampleName,...
    UM_limit_cosmo));

% Convert to polar stereographic for calculation measurement
[DV_limit_x, DV_limit_y] = ll2ps(DV_limit_lat, ...
    DV_limit_long);

[UM_limit_x, UM_limit_y] = ll2ps(UM_limit_lat, ...
    UM_limit_long);


%calculate euclidean distance to the Hatherton glacier centerline for each
%sample. This is probably a good enough approximation
for jj = 1:length(DV_limit_cosmo)
    DV_distance_to_centerline(jj) = min(sqrt((centerline_x_interp - ...
        DV_limit_x(jj)).^2 + (centerline_y_interp - DV_limit_y(jj)).^2));
    
end


for kk = 1:length(UM_limit_cosmo)
    UM_distance_to_centerline(kk) = min(sqrt((centerline_x_interp - ...
        UM_limit_x(kk)).^2 + (centerline_y_interp - UM_limit_y(kk)).^2));
   
end

%Concatenate all together, since we are fitting one curve to these two
%populations

DVBV_limit_distance_to_centerline = [DV_distance_to_centerline,...
    UM_distance_to_centerline];
DVBV_limit_elev = [DV_limit_elev', UM_limit_elev'];

%plot it up to make sure it looks reasonable
figure(3)
plot(DVBV_limit_distance_to_centerline, DVBV_limit_elev, 'ko'); hold on

% Good! Now fit a parabola to this, and evaluate it out to estimate the
% centerline thickness. We mirror elevations across the glacier to ensure
% a concave-down parabola, symmetric about the glacier centerline. Probably
% a little heavy-handed, but we are already making a biggish assumption
% that the cross surface profile is a parabola and that the glacier width
% does not change through time (ie., centerline stays put!).
n = 1; %2nd order polynomial
% DV_poly = polyfit([-DVBV_limit_distance_to_centerline,...
%     DVBV_limit_distance_to_centerline] , [DVBV_limit_elev, DVBV_limit_elev],n);

DV_poly = polyfit(DVBV_limit_distance_to_centerline, DVBV_limit_elev,n);
DV_surface_profile = polyval(DV_poly, linspace(0,...
    max(DVBV_limit_distance_to_centerline), 100));




% now plot up the surface profile to make sure it's a good fit
plot(linspace(0, max(DVBV_limit_distance_to_centerline), 100), DV_surface_profile)
title ('Dubris Valley')
xlabel ('Distance to centerline (m)')
ylabel ('Elevation (m)')
set (gca, 'fontsize', 13)

% assume the modern represents the largest difference between margin and
% centerline (100 m elevation in 2500 m distance) and the LGM projection
% you calculated here is the minimum difference (something like 131 m per
% 6800 m distance). Now for each sample, project it to the centerline using
% the average of these two.
data.output.distance_to_centerline = nan(size(data.output.t10St));

[sample_x, sample_y] = ll2ps(data.output.lat, data.output.lon);

for jj = 1:length(data.output.distance_to_centerline)
    data.output.distance_to_centerline(jj) = min(sqrt((centerline_x_interp - ...
        sample_x(jj)).^2 + (centerline_y_interp - sample_y(jj)).^2));
    
end

dS_dx = mean([100/2500, 131/6800]);
dS_dx = 100/2500;
dS_dx = 131/6824;

centerline_elev = data.output.elev + dS_dx.*data.output.distance_to_centerline; 
centerline_elev_err = dS_dx.*data.output.distance_to_centerline./2;

figure(5);%clf;  hold on; box on; grid on
DV_plot = errorbar(data.output.t10St(DV_index),centerline_elev(DV_index), centerline_elev_err(DV_index), 'ks');
BV_plot = errorbar(data.output.t10St(BV_index),centerline_elev(BV_index), centerline_elev_err(BV_index), 'ks');

DAN_plot = errorbar(data.output.t10St(DAN_index), centerline_elev(DAN_index), centerline_elev_err(DAN_index),'rs');
UM_plot = errorbar(data.output.t10St(UM_index), centerline_elev(UM_index), centerline_elev_err(UM_index),'rs');

plot(data.output.t10St(DV_index),centerline_elev(DV_index), 'ko')
plot(data.output.t10St(BV_index),centerline_elev(BV_index), 'ko')

plot(data.output.t10St(DAN_index), centerline_elev(DAN_index), 'rs')
plot(data.output.t10St(UM_index), centerline_elev(UM_index), 'rs', 'markerfacecolor', 'b')

xlim([0 12e3])
xlabel ('yr BP')
ylabel ('Projected glacier centerline elevation (m)')
set(gca, 'fontsize', 14)


%% Best estimate of LGM surface elevation is mean of highest erratics and
% linear fit at x = 0. The elevation span from the highest erratic to the
% mean could probably be considered 2 SDs, since the erratic is a strict
% minimum. Hence, the factor of 1.96 below in LGM_surface_*_sigma

LGM_surface_DV_elev = mean(polyval(DV_poly, [0,...
    DVBV_limit_distance_to_centerline(DVBV_limit_elev==max(DVBV_limit_elev))]));
LGM_surface_DV_sigma = 1/1.96 * abs(diff(polyval(DV_poly, [0,...
    DVBV_limit_distance_to_centerline(DVBV_limit_elev==max(DVBV_limit_elev))]))./2);

LGM_surface_MV_elev = mean(polyval(MV_poly, [0,...
    MV_limit_distance_to_centerline(MV_limit_elev==max(MV_limit_elev))]));
LGM_surface_MV_sigma = 1/1.96 * abs(diff(polyval(MV_poly, [0,...
    MV_limit_distance_to_centerline(MV_limit_elev==max(MV_limit_elev))]))./2);

LGM_surface_LW_elev = mean(polyval(LW_poly, [0,...
    LW_limit_distance_to_centerline(LW_limit_elev==max(LW_limit_elev))]));
LGM_surface_LW_sigma = 1/1.96 * abs(diff(polyval(LW_poly, [0,...
    LW_limit_distance_to_centerline(LW_limit_elev==max(LW_limit_elev))]))./2);

DVBV_shift = LGM_surface_DV_elev - mean(DV_limit_elev);
DAN_UM_shift = LGM_surface_DV_elev - mean(UM_limit_elev);
MV_walls_shift = LGM_surface_MV_elev - mean(MV_walls_limit_elev);
MV_floor_shift = LGM_surface_MV_elev - mean(MV_valley_limit_elev);
LW_valley_shift = LGM_surface_LW_elev - mean(LW_valley_limit_elev);
LW_walls_shift = LGM_surface_LW_elev - mean(LW_walls_limit_elev);


%% Convert to fraction of total elevation change, where 1 is LGM, 0 is modern glacier,
% then turn that into a centerline elevation
LW_margin = 870;
MV_margin = 1000;
DVBV_margin = 1130;

LW_x_P_ind = 13; 
MV_x_P_ind = 31; 
DAN_x_P_ind = 46; 

load('../../TEST_DH/output/DH_run_1a.mat', 'S_modern', 'S_modern2')

%LW
LW_algae_elev_frac = (LW_algae.Elevation - LW_margin)...
    ./(polyval(LW_poly, max(LW_limit_distance_to_centerline)) - LW_margin);

LW_walls_cosmo_elev_frac = (data.parsed_output.LW.erratics.elev(LW_walls_cosmo_index) ...
    - LW_margin) ./ (polyval(LW_poly,...
    LW_limit_distance_to_centerline(LW_limit_elev==max(LW_limit_elev))) ...
    - LW_margin);

LW_valley_cosmo_elev_frac = (data.parsed_output.LW.erratics.elev(~LW_walls_cosmo_index) ...
    - LW_margin) ./ (polyval(LW_poly, max(LW_limit_distance_to_centerline)) ...
    - LW_margin);

%turn into centerline elevation through time
LW_algae2center = LW_algae_elev_frac.*...
        (LGM_surface_LW_elev - S_modern2(LW_x_P_ind)) + ...
         S_modern2(LW_x_P_ind);
     
LW_walls_cosmo2center = LW_walls_cosmo_elev_frac.*...
        (LGM_surface_LW_elev - S_modern2(LW_x_P_ind)) + ...
         S_modern2(LW_x_P_ind);

LW_valley_cosmo2center = LW_valley_cosmo_elev_frac.*...
        (LGM_surface_LW_elev - S_modern2(LW_x_P_ind)) + ...
         S_modern2(LW_x_P_ind);

%MV elevation fraction
MV_walls_algae_elev_frac = (MV_walls_algae.Elevation - MV_margin)...
    ./(polyval(MV_poly, MV_limit_distance_to_centerline(MV_limit_elev == ...
    max(MV_limit_elev))) - MV_margin);

MV_floor_algae_elev_frac = (MV_floor_algae.Elevation - MV_margin)...
    ./(polyval(MV_poly, max(MV_limit_distance_to_centerline)) - MV_margin);

MV_walls_cosmo_elev_frac = (data.parsed_output.MV_walls.erratics.elev ...
    - MV_margin) ./ (polyval(MV_poly,...
    MV_limit_distance_to_centerline(MV_limit_elev == ...
    max(MV_limit_elev))) - MV_margin);

MV_floor_cosmo_elev_frac = (data.parsed_output.MV_floor.erratics.elev ...
    - MV_margin) ./ (polyval(MV_poly,...
    max(MV_limit_distance_to_centerline)) - MV_margin);

%MV centerline elevation through time
MV_walls_algae2center = MV_walls_algae_elev_frac.*...
        (LGM_surface_MV_elev - S_modern2(MV_x_P_ind)) + ...
         S_modern2(MV_x_P_ind);
MV_floor_algae2center = MV_floor_algae_elev_frac.*...
        (LGM_surface_MV_elev - S_modern2(MV_x_P_ind)) + ...
         S_modern2(MV_x_P_ind);
     
MV_walls_cosmo2center = MV_walls_cosmo_elev_frac.*...
        (LGM_surface_MV_elev - S_modern2(MV_x_P_ind)) + ...
         S_modern2(MV_x_P_ind);
MV_floor_cosmo2center = MV_floor_cosmo_elev_frac.*...
        (LGM_surface_MV_elev - S_modern2(MV_x_P_ind)) + ...
         S_modern2(MV_x_P_ind);


%DV, BV, DAN, UM
DV_cosmo_elev_frac = (data.parsed_output.DV.erratics.elev ...
    - DVBV_margin) ./ (polyval(DV_poly,...
    max(DVBV_limit_distance_to_centerline)) - DVBV_margin);
BV_cosmo_elev_frac = (data.parsed_output.BV.erratics.elev ...
    - DVBV_margin) ./ (polyval(DV_poly,...
    max(DVBV_limit_distance_to_centerline)) - DVBV_margin);

DAN_cosmo_elev_frac = (data.parsed_output.DAN.erratics.elev ...
    - DVBV_margin) ./ (polyval(DV_poly,...
    UM_distance_to_centerline(UM_limit_elev == ...
    max(UM_limit_elev))) - DVBV_margin);

UM_cosmo_elev_frac = (data.parsed_output.UM.erratics.elev ...
    - DVBV_margin) ./ (polyval(DV_poly,...
    UM_distance_to_centerline(UM_limit_elev == ...
    max(UM_limit_elev))) - DVBV_margin);

DVBV_algae_elev_frac = (DVBV_algae.Elevation - DVBV_margin) ./ ...
    (polyval(DV_poly, max(DVBV_limit_distance_to_centerline)) - DVBV_margin);

DAN_algae_elev_frac = (DAN_algae.Elevation - DVBV_margin) ./ ...
    (polyval(DV_poly,UM_distance_to_centerline(UM_limit_elev == ...
    max(UM_limit_elev))) - DVBV_margin);

%DV, BV, etc centerline elevation through time
DVBV_algae2center = DVBV_algae_elev_frac.*...
        (LGM_surface_DV_elev - S_modern2(DAN_x_P_ind)) + ...
         S_modern2(DAN_x_P_ind);
     
DAN_algae2center = DAN_algae_elev_frac.*...
        (LGM_surface_DV_elev - S_modern2(DAN_x_P_ind)) + ...
         S_modern2(DAN_x_P_ind);
     
DV_cosmo2center = DV_cosmo_elev_frac.*...
        (LGM_surface_DV_elev - S_modern2(DAN_x_P_ind)) + ...
         S_modern2(DAN_x_P_ind);
     
BV_cosmo2center = BV_cosmo_elev_frac.*...
        (LGM_surface_DV_elev - S_modern2(DAN_x_P_ind)) + ...
         S_modern2(DAN_x_P_ind);
     
DAN_cosmo2center = DAN_cosmo_elev_frac.*...
        (LGM_surface_DV_elev - S_modern2(DAN_x_P_ind)) + ...
         S_modern2(DAN_x_P_ind);
     
UM_cosmo2center = UM_cosmo_elev_frac.*...
        (LGM_surface_DV_elev - S_modern2(DAN_x_P_ind)) + ...
         S_modern2(DAN_x_P_ind);

% %%
 save geochron_constraints.mat LGM_surface*_elev LGM_surface_*_sigma ...
     *_shift *elev_frac LW_walls_cosmo_index *cosmo2center *algae2center ...
     S_modern S_modern2 *x_P_ind -append