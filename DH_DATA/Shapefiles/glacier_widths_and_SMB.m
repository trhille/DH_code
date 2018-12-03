function [width, SMB_mean, distance_along_centerline] = glacier_widths_and_SMB(RACMO_paths, flowline_path, buffer_path, nodes)
%glacier_widths_and_SMB automatically measures widths of a polygon (glacier outline
%shapefile) perpendicular to a flowline (also a shapefile), and calculates
%the width-averaged surface mass balance at each node
%  
%RACMO_paths is a cell array with paths to precip, subl, and runoff IN THAT
%ORDER
%flowline = Polar stereographic pairs of flowline coordinates, starting at
%the glacier mouth. Use Qgis to make this, then use interp1 to densify.
%Read into matlab with shaperead
%buffer = the polygon of the glacier outline, made in qgis. use shaperead
%nodes = number of evenly-spaced points at which to measure width.
%% get shapefiles into usable format. Buffer is fine the way it is.
addpath '/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton/Data/RACMO/RACMO2.3';
addpath '/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton/Data/RACMO/RACMO2.1_ADL055';
addpath('/Users/trevorhillebrand/Documents/MATLAB/Toolboxes/AntarcticMappingTools_v5.00/AntarcticMappingTools');

if iscell(RACMO_paths)  %if you have three different .nc files with precip, sublimation, and runoff
    
    if strfind(RACMO_paths{1}, '.nc')
    precip = RACMO_nc2mat(RACMO_paths{1});
    subl = RACMO_nc2mat(RACMO_paths{2});
    runoff = RACMO_nc2mat(RACMO_paths{3});
    
    elseif strfind(RACMO_paths{1}, '.mat')
    load(RACMO_paths{1});
    load(RACMO_paths{2});
    load(RACMO_paths{3});
    SMB = precip.precip + subl.subl - runoff.runoff; 
    [SMB_x, SMB_y] = ll2ps(precip.lat, precip.lon);
    
elseif ~iscell(RACMO_paths) %if you have just one .nc file with just total yearly SMB
    
    RACMO_SMB = RACMO_nc2mat(RACMO_paths);
    SMB = mean(RACMO_SMB.smb, 4); %RACMO_nc2mat loads the smb as a 4-D array
    [SMB_x, SMB_y] = ll2ps(RACMO_SMB.lat, RACMO_SMB.lon);
   
end

buffer = shaperead(buffer_path);

flowline_struct = shaperead(flowline_path);
flowline_struct.X = flowline_struct.X(1:end-1);
flowline_struct.Y = flowline_struct.Y(1:end-1);

flowx = linspace(flowline_struct.X(1), flowline_struct.X(end), nodes);
flowy = interp1(flowline_struct.X, flowline_struct.Y, flowx);
flowline = [flowx; flowy]';

%% Normal lies in the null space of the matrix A - B, where A and B are each
%(x,y) pairs
%calculate normals to left and right of start and end of flowline (looking
%upstream)
temp_halfwidth = 5e4; %this just needs to be wider than any point along the flowline. 50 km is a good halfwidth

Lnormal_start = flowline(1,:) + temp_halfwidth.*(null(flowline(2,:) - flowline(1,:)))';
Rnormal_start = flowline(1,:) - temp_halfwidth.*(null(flowline(2,:) - flowline(1,:)))';
Lnormal_end = flowline(end,:) + temp_halfwidth.*(null(flowline(end,:) - flowline(end-1,:)))';
Rnormal_end = flowline(end,:) - temp_halfwidth.*(null(flowline(end,:) - flowline(end-1,:)))';

Lnormal = Lnormal_start;
Rnormal = Rnormal_start;

for jj = 2:size(flowline,1)-1
Lnormal_temp = flowline(jj, :) + temp_halfwidth.*(null(flowline(jj+1,:)-flowline(jj-1,:)))';
Rnormal_temp = flowline(jj, :) - temp_halfwidth.*(null(flowline(jj+1,:)-flowline(jj-1,:)))';

Lnormal = [Lnormal; Lnormal_temp];
Rnormal = [Rnormal; Rnormal_temp];
end

Lnormal = [Lnormal; Lnormal_end];
Rnormal = [Rnormal; Rnormal_end];

%% Now find the points where these normal vectors intersect the buffer
count = 1;
SMB_mean = zeros(size(flowline,1), 1);
for jj = 1:size(flowline,1)
  
    normalx = [Rnormal(jj,1), Lnormal(jj,1)];
    normaly = [Rnormal(jj,2), Lnormal(jj,2)];
    
    % This gives the x and y coordinates of the points where widthlines
    % intersect the buffer
    [xx,yy] = polyxpoly(normalx, normaly, buffer.X, buffer.Y); 
    
    plot(xx,yy)
    hold on
    
    %The width is then the euclidian distance between these two points
   
    if isempty(xx) == 0 && isempty(yy) == 0 
        width(count) = pdist2([xx(1) yy(1)],[xx(2) yy(2)]);
   
    
    %temporarily fill in between these two points to calculate
    %width-averaged SMB

    xx_interp = linspace(xx(1), xx(2), ceil(width(count)/1000));
    yy_interp = interp1(xx, yy, xx_interp);
    
    SMB_profile = griddata(SMB_x, SMB_y, SMB, xx_interp, yy_interp, 'linear')./917; %interpolate and convert to m ice equivalent (specific gravity of ice =0.917 and convert mm to m)
    SMB_mean(jj) = sum(SMB_profile)./length(SMB_profile);
  
    
    else
        width(count) = NaN;
    end
    count = count+1;
end

plot(buffer.X, buffer.Y)
plot(flowline(:,1), flowline(:,2), 'k', 'linewidth', 2)
daspect([1 1 1])

distance_vec = [0];
for jj = 2:size(flowline,1)
distance_temp= [pdist2([flowline(jj-1,1) flowline(jj-1,2)], [flowline(jj,1) flowline(jj,2)])'];

distance_vec  = [distance_vec; distance_temp];
end

distance_along_centerline = cumsum(distance_vec);
keyboard
end






