function [ b_dot_edges, b_dot_P ] = ...
    iterate_bdot( glacier_surface, precip_at_sl, lapse, x_P, x_edges)
%iterate_bdot is used to update the surface mass balance at
%   Detailed explanation goes here


  b_dot_P = precip_at_sl + lapse.*glacier_surface; 
  % keep dependence on surface negative, following Patankar (1980) pg 48
    
 
% interpolate onto edges of mesh
% ------------------------------ 
   b_dot_edges = interp1( x_P, b_dot_P, x_edges );  

end

