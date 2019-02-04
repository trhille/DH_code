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
n = 2; %2nd order polynomial
LW_poly = polyfit([-LW_limit_distance_to_centerline,LW_limit_distance_to_centerline] , [LW_limit_elev, LW_limit_elev],n);

LW_surface_profile = polyval(LW_poly, linspace(0, max(LW_limit_distance_to_centerline), 100));

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
n = 2; %2nd order polynomial
MV_poly = polyfit([-MV_limit_distance_to_centerline,...
    MV_limit_distance_to_centerline] , [MV_limit_elev, MV_limit_elev],n);

MV_surface_profile = polyval(MV_poly, linspace(0,...
    max(MV_limit_distance_to_centerline), 100));

% now plot up the surface profile to make sure it's a good fit
plot(linspace(0, max(MV_limit_distance_to_centerline), 100), MV_surface_profile)
title ('Magnis Valley')
xlabel ('Distance to centerline (m)')
ylabel ('Elevation (m)')
set (gca, 'fontsize', 13)

%% Dubris Valley and Updog Mountain
%Also copy-pasted from above, then edited. 

%% Lake Wellman

%define which samples consitute the LGM limit in the valley and on the
%walls nearest the glacier
LW_valley_limit_cosmo = {'13-HAT-030-LW', '13-HAT-029-LW', '13-HAT-031-LW'};
LW_walls_limit_cosmo = {'13-HAT-044-LW', '13-HAT-047-LW', '13-HAT-006-LW',...
    '13-HAT-007-LW'};
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
n = 2; %2nd order polynomial
LW_poly = polyfit([-LW_limit_distance_to_centerline,LW_limit_distance_to_centerline] , [LW_limit_elev, LW_limit_elev],n);

LW_surface_profile = polyval(LW_poly, linspace(0, max(LW_limit_distance_to_centerline), 100));

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
n = 2; %2nd order polynomial
DV_poly = polyfit([-DVBV_limit_distance_to_centerline,...
    DVBV_limit_distance_to_centerline] , [DVBV_limit_elev, DVBV_limit_elev],n);

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


centerline_elev = data.output.elev + mean([100/2500, 131/6800]).*data.output.distance_to_centerline; 
centerline_elev_err = mean([100/2500, 131/6800]).*data.output.distance_to_centerline./2;

figure(5);clf;  hold on; box on; grid on
DV_plot = errorbar(data.output.t10St(DV_index),centerline_elev(DV_index), centerline_elev_err(DV_index), 'ks');
BV_plot = errorbar(data.output.t10St(BV_index),centerline_elev(BV_index), centerline_elev_err(BV_index), 'ks');

DAN_plot = errorbar(data.output.t10St(DAN_index), centerline_elev(DAN_index), centerline_elev_err(DAN_index),'rs');
UM_plot = errorbar(data.output.t10St(UM_index), centerline_elev(UM_index), centerline_elev_err(UM_index),'rs');

% plot(data.output.t10St(DV_index),centerline_elev(DV_index), 'ko')
% plot(data.output.t10St(BV_index),centerline_elev(BV_index), 'ko')
% 
% plot(data.output.t10St(DAN_index), centerline_elev(DAN_index), 'rs')
% plot(data.output.t10St(UM_index), centerline_elev(UM_index), 'rs')

xlim([0 12e3])
xlabel ('yr BP')
ylabel ('Projected glacier centerline elevation (m)')
set(gca, 'fontsize', 14)
legend([DV_plot, DAN_plot], {'Dubris and Bibra Valleys', 'Danum Platform and Updog Mtn'})



