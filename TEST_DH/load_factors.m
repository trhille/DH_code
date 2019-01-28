function [ E_P, E_w, E_e, ...
           fs_P, fs_w, fs_e, ...
           scaling_P, scaling_w, scaling_e] = load_factors( x_P, x_w, x_e, dx_P, dx_w, dx_e, S_modern, B_P, W_P )


%--------------------------------------------------------------
%
%  set up steady sliding values slip_P(x) at positions x_P
%  where surface S(x,t) is defined
%  interpolate bed to interval midpoints where fluxes are
%  calculated
%
%
%  Following equation in the form (Budd 1979; Oerlemans 2001; 
%                                  Cuffey and Paterson 2010, pg. 295+, 339+,
%                                  464+ )
%      ubar = ubar_deformation + ubar_sliding
%      shear stress tau_d = rho * g * H * alpha  
%      assume tau_d = tau_b
%      ubar = (fd * H * tau_d^n) + (fs * (1/H) * tau_d^m)
%      fd = 2 * E * A0 * exp(-Q/R*T) in Pa^-3 yr^-1; prescribe E
%      fs = prescribe factor in Pa^3 m^2 yr^-1
%      assume n = m = 3!!
%
%---------------------------------------------------------------


   
global s_per_year


% Deformation enhancement factor:
% -------------------------------
   E_P    =  1 * ones( size(x_P) );  % Set equal to one there is no enhancement
                                     % Scaling factor that has no units
                                     
%    index_50km = find(x_P<=50000, 1, 'last');                                  
%    E_P(1:index_50km) = 0;
                                 
                     
% % Values from min search
% % ======================
% disp('Using best min search factor for E!')
% load min7.mat E_P  
% disp('loading E_P from min7.mat')
% fd_P = interp1(x_P_best, factor_P_best, x_P, 'linear', 'extrap');
 
  [E_w, E_e ] = get_edge_values_quadratic ( E_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );
  
   

% Sliding factor:
% ---------------
% Make sure in units Pa^3 m^2 yr^-1

  % fs_P = (5.7e-20 * s_per_year) * ones(size(x_P));  % Budd et al. (1979) but check other values?
 
    fs_P = 0.5e-11 * ones(size(x_P)); 
   
  %  fs_P(1:index_50km) = 1e-11;
    

  [fs_w, fs_e ] = get_edge_values_quadratic ( fs_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );


  
  
  % Scaling factor based on geometry, used to scale ice flux
  % ---------------------------------------------------------
  load DH_DATA/DH_cross_sect_area.mat Darwin_cross_sect_area...
      Darwin_distance_along_centerline
  
  load DH_accum_width_velocity.mat Darwin_width_x Darwin_width_values
  
  W_P_true = interp1(Darwin_width_x, Darwin_width_values, x_P, 'linear', 'extrap');
  
  cross_sect_area = interp1(Darwin_distance_along_centerline + x_P(1),...
      Darwin_cross_sect_area, x_P, 'linear', 'extrap');
  
  scaling_P = cross_sect_area./((S_modern - B_P).*W_P_true);
  
  [ scaling_w, scaling_e ] = get_edge_values_quadratic (scaling_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );

  
                           