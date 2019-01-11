%calculate tributary fluxes at Darwin and Hatherton Glaciers using flux
%gates, ice thickness from Mette, and MEASUREs velocities.
Darwin_data_path = ['/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton', ...
    '/Modeling/Koutnik model/DH_code/DH_DATA/Shapefiles/Darwin Flowline/'];

Hatherton_data_path = ['/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton', ...
    '/Modeling/Koutnik model/DH_code/DH_DATA/Shapefiles/Hatherton Flowline/'];

Darwin_1 = [Darwin_data_path, 'Darwin_tributary_1_fluxgate.xlsx'];
Darwin_2 = [Darwin_data_path, 'Darwin_tributary_2_fluxgate.xlsx'];
Darwin_3 = [Darwin_data_path, 'Darwin_tributary_3_fluxgate.xlsx'];

Hatherton_1 = [Hatherton_data_path, 'Hatherton_tributary1_fluxgate.xlsx'];
Hatherton_2 = [Hatherton_data_path, 'Hatherton_tributary2_fluxgate.xlsx'];


%% load data for fluxgates
Darwin_icethick_1 = readtable(Darwin_1, 'Sheet', 'Ice_Thickness');
Darwin_vel_1 = readtable(Darwin_1, 'Sheet', 'FlowVelocity');
Darwin_icethick_2 = readtable(Darwin_2, 'Sheet', 'Ice_Thickness');
Darwin_vel_2 = readtable(Darwin_2, 'Sheet', 'FlowVelocity');
Darwin_icethick_3 = readtable(Darwin_3, 'Sheet', 'Ice_Thickness');
Darwin_vel_3 = readtable(Darwin_3, 'Sheet', 'FlowVelocity');

Hatherton_icethick_1 = readtable(Hatherton_1, 'Sheet', 'Ice_Thickness');
Hatherton_vel_1 = readtable(Hatherton_1, 'Sheet', 'FlowVelocity');
Hatherton_icethick_2 = readtable(Hatherton_2, 'Sheet', 'Ice_Thickness');
Hatherton_vel_2 = readtable(Hatherton_2, 'Sheet', 'FlowVelocity');

%% Calculate flux

%calculate cross sectional area for each flux gate
Darwin_1_cross_sectional_area = gradient(Darwin_icethick_1.Distance_m).* ...
    Darwin_icethick_1.Ice_Thickness_m;
Darwin_2_cross_sectional_area = gradient(Darwin_icethick_2.Distance_m).* ...
    Darwin_icethick_2.Ice_Thickness_m;
Darwin_3_cross_sectional_area = gradient(Darwin_icethick_3.Distance_m).* ...
    Darwin_icethick_3.Ice_Thickness_m;

Hatherton_1_cross_sectional_area = gradient(Hatherton_icethick_1.Distance_m).* ...
    Hatherton_icethick_1.Ice_Thickness_m;
Hatherton_2_cross_sectional_area = gradient(Hatherton_icethick_2.Distance_m).* ...
    Hatherton_icethick_2.Ice_Thickness_m;

%interpolate velocities onto thickness grid
Darwin_1_velocity_interp = interp1(Darwin_vel_1.Distance_m,...
    Darwin_vel_1.FlowVelocity_m_yr, Darwin_icethick_1.Distance_m);
Darwin_2_velocity_interp = interp1(Darwin_vel_2.Distance_m,...
    Darwin_vel_2.FlowVelocity_m_yr, Darwin_icethick_2.Distance_m);
Darwin_3_velocity_interp = interp1(Darwin_vel_3.Distance_m,...
    Darwin_vel_3.FlowVelocity_m_yr, Darwin_icethick_3.Distance_m);

Hatherton_1_velocity_interp = interp1(Hatherton_vel_1.Distance_m,...
    Hatherton_vel_1.FlowVelocity_m_yr, Hatherton_icethick_1.Distance_m);
Hatherton_2_velocity_interp = interp1(Hatherton_vel_2.Distance_m,...
    Hatherton_vel_2.FlowVelocity_m_yr, Hatherton_icethick_2.Distance_m);

%so calculate that flux
Darwin_1_flux = sum(Darwin_1_cross_sectional_area.*Darwin_1_velocity_interp);
Darwin_2_flux = sum(Darwin_2_cross_sectional_area.*Darwin_2_velocity_interp);
Darwin_3_flux = sum(Darwin_3_cross_sectional_area.*Darwin_3_velocity_interp);

Hatherton_1_flux = sum(Hatherton_1_cross_sectional_area.*Hatherton_1_velocity_interp);
Hatherton_2_flux = sum(Hatherton_2_cross_sectional_area.*Hatherton_2_velocity_interp);
    