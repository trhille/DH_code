function [ output_args ] = untitled( input_args )
%UNTITLED Summary of this function goes here
%   This function estimates tributary flux for Darwin and Hatherton
%   Glaciers. It needs a RACMO .nc file (or .mat file made by
%   glacier_widths_and_SMB.m) and a shapefile for each tributary.
%% get shapefiles into usable format. Buffer is fine the way it is.
addpath '/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton/Data/RACMO/RACMO2.3';
addpath '/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton/Data/RACMO/RACMO2.1_ADL055';
addpath('/Users/trevorhillebrand/Documents/MATLAB/Toolboxes/AntarcticMappingTools_v5.00/AntarcticMappingTools');

Hat_trib1_path = '/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton/Modeling/Koutnik model/DH_code/DH_DATA/Shapefiles/Hatherton Flowline/Hatherton_tributary1.shp';
Hat_trib2_path = '/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton/Modeling/Koutnik model/DH_code/DH_DATA/Shapefiles/Hatherton Flowline/Hatherton_tributary2.shp';

Hat_trib1 = shaperead(Hat_trib1_path);
Hat_trib2 = shaperead(Hat_trib2_path);

%% Now load RACMO data and reproject it onto a 1 km grid. Do you need to divide that values by the resolution, or are they per km^2?

load precip.mat
load subl.mat
end

