% run_min_search_E

E_vec    = [ 0.1 0.5 1 5 10 15 20 30 ];
% E_vec = [1]; % 3 10];
% fs_vec = [ 0.1e-11 0.2e-11 0.3e-11 0.4e-11 0.5e-11 0.6e-11 0.7e-11 0.8e-11 0.9e-11 1e-11];
fs_vec = [ 1e-12 1e-11 1e-10 ];


N_E           = length( E_vec );
N_fs          = length( fs_vec );


% Still find best parameter set for Darwin, then for Hatherton
% But now search over combinations of E and sliding, include E=1 and fs=0?

% CHANGE iterations = 1 or = 2 here!!
% Might be working to run through both together but don't have file saving
% setup and don't have it double checked. 

% for iterations = 1   % Minimize for both Darwin (1) and Hatherton (2)
 for iterations = 2    
    
    if (iterations == 1)
        disp('Running for Darwin...')
        length_run      = length(x_P);
        save_mismatch   = NaN * ones(length_run, N_E*N_fs);
        save_E          = NaN * ones(1, length_run);
        save_fs         = NaN * ones(1, length_run);
        E_P_orig        = E_P;
        fs_P_orig       = fs_P;
        Darwin_velocity = interp1(Darwin_measures_centerline_distance, Darwin_measures_flowspeed, [x_w(1) x_e]);
        Darwin_surface  = S_modern;
        RMS_mismatch_matrix = NaN * ones(N_E, N_fs);
        
    

 for xpos = 2:length_run-1   % loop over all x positions -- upstream
    
disp(['Evaluating E and fs for x-position ',int2str(xpos), ' of total=',int2str(length_run)]);

RMS_mismatch_matrix = NaN * ones( N_E, N_fs );


%  Loop over E_x
 for i_E = 1:N_E
     
      E_use = E_vec(i_E);
      E_P(xpos) = E_use;
      
    % recalculate values on the edges.
      [ E_w, E_e ] = get_edge_values_quadratic ...
                               ( E_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );
                       
                           
   for i_fs = 1:N_fs                       
         
           
    fs_use = fs_vec(i_fs);

    fs_P(xpos) = fs_use;
       
    % recalculate values on the edges.
      [fs_w, fs_e ] = get_edge_values_quadratic ...
                               ( fs_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );
    
  % Still get some negatives!                          
  fs_w(find(fs_w < 0)) = 1e-11;                         
  fs_e(find(fs_e < 0)) = 1e-11;                         
  E_w(find(E_w < 0))   = 1;
  E_e(find(E_e < 0))   = 1;
  
                           
  [ S_P(time,:), h_P(time,:), ...  
    dS_dx_P_xt(time,:),  ...
    dS_dx_edges_xt(time,:) ] = calc_h_0( x_P, x_w, x_e, dx_P, B_P, B_w, B_e, ...
                                         W_P, W_w, W_e, b_dot_edges(time,:), ...
                                         E_w, E_e, fs_w, fs_e, ...
                                         Q_0_in, S_at_GL(time), ...
                                         A_eff_edges_xt(time,1:end-1), ...
                                         A_eff_edges_xt(time,2:end), ...
                                         flux_add_w, flux_add_e, ...
                                         deformation_only, deformation_plus_sliding, sliding_only);
                                                         
% edge values of ice thickness
% ============================
  [ h_w(time,:), ...
    h_e(time,:) ] = get_edge_values_quadratic( h_P(time,:), x_P, x_w, x_e, ...
                                               dx_P, dx_w, dx_e );
              
% Find the flux
% ==============
 flux_kin_P_xt(time,:) = calc_flux_kin( x_P, x_P, dx_P, W_P, ...
                                        h_dot(time,:), b_dot_P(time,:), Q_0_in );
                                                                                             
[ flux_w, flux_e ] = get_edge_values_quadratic ( flux_kin_P_xt(time,:), ...
                                                 x_P, x_w, x_e, ...
                                                 dx_P, dx_w, dx_e );
   
 flux_edges_kin_xt(time,:) = [ Q_0_in flux_e];   % Need to replace with Q_0_in, not extrapolated value.
                   
 
 flux_edges_dyn_xt(time,:) = calc_flux_dyn( x_edges, [h_w(time,1) h_e(time,:)], ...
                                            dS_dx_edges_xt(time,:), ...
                                            [E_w(1) E_e], [fs_w(1) fs_e], ...
                                            [W_w(1) W_e], A_eff_edges_xt(time,:), ...
                                            deformation_only, deformation_plus_sliding, sliding_only);   
 flux_edges_dyn_xt(1,1)    = Q_0_in; 
           
 
 % Velocity
 % ========
 surf_vel_estimate = (5/4) * abs(flux_edges_kin_xt(1,:)) ./ ([W_w(1) W_e] .* [h_w(1,1) h_e(1,:)]);
  
 
 % Residual
 % ========
    std_dev = 1; 

% %  % Surface profile:
%   residual = sqrt( mean( ( (abs(Darwin_surface) - abs(S_P(1,:)))/std_dev ).^2 ) );           
 %   residual = sqrt( mean( ( (abs(Darwin_surface(xpos)) - abs(S_P(1,xpos)))/std_dev ).^2 ) );  

 % % Surface velocity:
  % residual = sqrt( mean( ( (abs(Darwin_velocity) - abs(surf_vel_estimate))/std_dev ).^2 ) );  
  %  residual = abs(Darwin_velocity(xpos) - surf_vel_estimate(xpos));  % local mismatch?
 
 % % Both:
 %  residual = sqrt( mean( ( (abs(Darwin_velocity) - abs(surf_vel_estimate))/std_dev ).^2 ) ) + ...
 %             sqrt( mean( ( (abs(Darwin_surface) - abs(S_P(1,:)))/std_dev ).^2 ) );    
   residual = sqrt( mean( ( (abs(Darwin_velocity(xpos)) - abs(surf_vel_estimate(xpos)))/std_dev ).^2 ) ) + ...
              sqrt( mean( ( (abs(Darwin_surface(xpos)) - abs(S_P(1,xpos)))/std_dev ).^2 ) );   
          
  
if (isreal(S_P(1,:)) == 0)
    residual = NaN;   % Don't pick from this value if gives imaginary.
end
  
   RMS_mismatch_matrix(i_E,i_fs) = residual;
   
   
   end % loop over fs values
   
  
  end  % loop over E values
    
    
   RMS_mismatch_matrix;
   min_value = min(min(RMS_mismatch_matrix));
  [minRow, minCol] = find(RMS_mismatch_matrix == min_value);
 
   
%    if (length(index)>1)  % Multiple values of E give same mismatch
%        disp('Multiple values of fs give same mismatch value -- set fs=1e-11')
%        index = find(fs_vec == 1e-11);
%    end


   RMS_best               = RMS_mismatch_matrix(minRow, minCol);
   best_fs                = fs_vec(minCol);
   best_E                 = E_vec(minRow);
   save_fs(xpos)          = best_fs;   
   save_E(xpos)           = best_E;
   save_mismatch(xpos, :) = RMS_mismatch_matrix(minRow, minCol);

   
 % Keep best value:
   fs_P(xpos) = best_fs; 
                                  
   E_P(xpos) = best_E;
 
     
 
end  % loop over xpositions

      x_P_min = x_P;
      x_w_min = x_w;
      x_e_min = x_e;
      
     % Get around issue that interpolating to edge points can give negative 
      save_E(1)   = save_E(2);
      save_E(end) = save_E(end-1);
      E_P_min     = save_E;
      
    % recalculate values on the edges.
      [ E_w_min, E_e_min ] = get_edge_values_quadratic ...
                               ( E_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );

     save_fs(1)   = save_fs(2);
     save_fs(end) = save_fs(end-1);
     fs_P_min     = save_fs;

    % recalculate values on the edges.
      [ fs_w_min, fs_e_min ] = get_edge_values_quadratic ...
                               ( fs_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );
  
 % Turns out that negatives can still come in. Just reset.                          
   fs_w_min(1)   = fs_w_min(2);
   fs_e_min(end) = fs_e_min(end-1);
   E_w_min(1)    = E_w_min(2);
   E_e_min(end)  = E_e_min(end-1);
   
   save min_Darwin_E_and_fs.mat E_P_min E_w_min E_e_min ...
                             fs_P_min fs_w_min fs_e_min ...
                             x_P_min x_w_min x_e_min  
                        
  
 [ S_P(time,:), h_P(time,:), ...  
    dS_dx_P_xt(time,:),  ...
    dS_dx_edges_xt(time,:) ] = calc_h_0( x_P, x_w, x_e, dx_P, B_P, B_w, B_e, ...
                                         W_P, W_w, W_e, b_dot_edges(time,:), ...
                                         E_w_min, E_e_min, fs_w_min, fs_e_min, ...
                                         Q_0_in, S_at_GL(time), ...
                                         A_eff_edges_xt(time,1:end-1), ...
                                         A_eff_edges_xt(time,2:end), ...
                                         flux_add_w, flux_add_e, ...
                                         deformation_only, deformation_plus_sliding, sliding_only);
                                                         
              
% Find the flux
% ==============
 flux_kin_P_xt(time,:) = calc_flux_kin( x_P, x_P, dx_P, W_P, ...
                                        h_dot(time,:), b_dot_P(time,:), Q_0_in );
                                                                                             
[ flux_w, flux_e ] = get_edge_values_quadratic ( flux_kin_P_xt(time,:), ...
                                                 x_P, x_w, x_e, ...
                                                 dx_P, dx_w, dx_e );
   
 flux_edges_kin_xt(time,:) = [ Q_0_in flux_e];   % Need to replace with Q_0_in, not extrapolated value.
                   
 
 % Velocity
 % ========
 surf_vel_estimate = (5/4) * abs(flux_edges_kin_xt(1,:)) ./ ([W_w(1) W_e] .* [h_w(1,1) h_e(1,:)]);
  
                          
 

figure(100)   % SURFACE and BED
subplot(2,1,1), plot(x_P/1000, B_P,'r')
hold on
plot(x_P/1000, S_modern, 'linewidth', 2)
plot(x_P/1000, S_P(1,:),'c--')
legend('Bed', 'Measured surface', 'Calculated surface', 'location', 'northwest')
title('Surface and Bed topography: DARWIN')
%xlabel('Distance along flowband (km)')
ylabel('Elevation (m)')
xlim([x_P(1) x_P(end)]/1000)
subplot(2,1,2)
plot(x_P/1000, S_modern-S_P(1,:), 'linewidth', 2)
hold on
plot([x_P(1) x_P(end)]/1000, [0 0],'k--')
title('Modern surface - calculated surface')
xlabel('Distance along flowband (km)')
ylabel('Elevation difference (m)')
xlim([x_P(1) x_P(end)]/1000)                            
      
      
figure(300)   % SURFACE VELOCITY 
subplot(2,1,1), plot(x_edges/1000, surf_vel_estimate, 'k--', 'linewidth', 2)
hold on
plot(Darwin_measures_centerline_distance/1000, Darwin_measures_flowspeed,'m', 'linewidth', 2)
title('Surface velocity: DARWIN')
ylabel('Velocity (m/yr)')
legend('Calculated velocity',  'MEaSUREs velocity', 'location', 'best')

measures_on_xedges = interp1( Darwin_measures_centerline_distance, Darwin_measures_flowspeed, x_edges, 'linear', 'extrap');

subplot(2,1,2), plot(x_edges/1000, measures_on_xedges - surf_vel_estimate*(5/4), 'b', 'linewidth', 2)
hold on
plot([x_edges(1) x_edges(end)]/1000, [0 0],'k--')
title('Observed - calculated surface velocity')
xlabel('Distance along flowband (km)')
ylabel('Velocity difference (m/yr)')

   
   




 % ------------------------------------------------------------------------                          
 % ------------------------------------------------------------------------                          
  elseif (iterations == 2)
        disp('Running for Hatherton...')
        length_run2      = length(x_P2);
        save_mismatch2   = NaN * ones(length_run2, N_E);
        save_E2          = NaN * ones(1, length_run2);
        save_fs2         = NaN * ones(1, length_run2);
        E_P_orig2        = E_P2;
        Hat_velocity    = interp1( Hat_measures_centerline_distance, Hat_measures_flowspeed, x_edges2, 'linear', 'extrap');
        Hat_surface     = S_modern2;
        RMS_mismatch_matrix2 = NaN * ones(N_E, N_fs);
        
    

 for xpos2 = 2:length_run2-1   % loop over all x positions -- upstream; 
                               % ignore edge values in minimization so don't get negative in interpolation
    
disp(['Evaluating E and fs for x-position ',int2str(xpos2), ' of total=',int2str(length_run2)]);

RMS_mismatch_matrix2 = NaN * ones( N_E, N_fs );


%  Loop over E_x
 for i_E = 1:N_E
     
      E_use = E_vec(i_E);
      E_P2(xpos2) = E_use;
      
    % recalculate values on the edges.
      [ E_w2, E_e2 ] = get_edge_values_quadratic ...
                               ( E_P2, x_P2, x_w2, x_e2, dx_P2, dx_w2, dx_e2 );
                       
                           
   for i_fs = 1:N_fs                       
         
           
    fs_use2 = fs_vec(i_fs);

    fs_P2(xpos2) = fs_use2;
       
    % recalculate values on the edges.
      [fs_w2, fs_e2 ] = get_edge_values_quadratic ...
                               ( fs_P2, x_P2, x_w2, x_e2, dx_P2, dx_w2, dx_e2 );
    
                           
 % Still possible to get some negatives! Not a perfect fix, but something to get past edge effect.                
  fs_w2(find(fs_w2 < 0)) = 1e-11;                         
  fs_e2(find(fs_e2 < 0)) = 1e-11;                         
  E_w2(find(E_w2 < 0))   = 1;
  E_e2(find(E_e2 < 0))   = 1;
  
                         
  [ S_P2(time2,:), h_P2(time,:), ...  
    dS_dx_P_xt2(time2,:),  ...
    dS_dx_edges_xt2(time2,:) ] = calc_h_0( x_P2, x_w2, x_e2, dx_P2, B_P2, B_w2, B_e2, ...
                                         W_P2, W_w2, W_e2, b_dot_edges2(time2,:), ...
                                         E_w2, E_e2, fs_w2, fs_e2, ...
                                         Q_0_in2, S_at_GL2(time2), ...
                                         A_eff_edges_xt2(time2,1:end-1), ...
                                         A_eff_edges_xt2(time2,2:end), ...
                                         flux_add_w2, flux_add_e2, ...
                                         deformation_only2, deformation_plus_sliding2, sliding_only2);
    
                                     
% edge values of ice thickness
% ============================
  [ h_w2(time2,:), ...
    h_e2(time2,:) ] = get_edge_values_quadratic( h_P2(time2,:), x_P2, x_w2, x_e2, ...
                                               dx_P2, dx_w2, dx_e2 );
              
% Find the flux
% ==============
 flux_kin_P_xt2(time2,:) = calc_flux_kin( x_P2, x_P2, dx_P2, W_P2, ...
                                        h_dot2(time2,:), b_dot_P2(time2,:), Q_0_in2 );
                                                                                             
[ flux_w2, flux_e2 ] = get_edge_values_quadratic ( flux_kin_P_xt2(time2,:), ...
                                                 x_P2, x_w2, x_e2, ...
                                                 dx_P2, dx_w2, dx_e2 );
   
 flux_edges_kin_xt2(time2,:) = [ Q_0_in2 flux_e2];   % Need to replace with Q_0_in, not extrapolated value.
                            
 
 
 flux_edges_dyn_xt2(time2,:) = calc_flux_dyn( x_edges2, [h_w2(time2,1) h_e2(time2,:)], ...
                                            dS_dx_edges_xt2(time2,:), ...
                                            [E_w2(1) E_e2], [fs_w2(1) fs_e2], ...
                                            [W_w2(1) W_e2], A_eff_edges_xt2(time2,:), ...
                                            deformation_only2, deformation_plus_sliding2, sliding_only2);           
                                     
 flux_edges_dyn_xt2(1,1)    = Q_0_in; 
                                        
         
 
 % Velocity
 % ========
 surf_vel_estimate2 = (5/4) * abs(flux_edges_kin_xt2(1,:)) ./ ([W_w2(1) W_e2] .* [h_w2(1,1) h_e2(1,:)]);
 average_vel_estimate2 = abs(flux_edges_kin_xt2(1,:)) ./ ([W_w2(1) W_e2] .* [h_w2(1,1) h_e2(1,:)]);
 
 
 % Residual
 % ========
    std_dev = 1; 

 % % Surface profile:
 %  residual2 = sqrt( mean( ( (abs(Hat_surface) - abs(S_P2(1,:)))/std_dev ).^2 ) );           

 % Surface velocity:
  % residual2 = sqrt( mean( ( (abs(Hat_velocity) - abs(surf_vel_estimate2))/std_dev ).^2 ) );  
  %  residual2 = abs(Hat_velocity(xpos2) - surf_vel_estimate2(xpos2));  % local mismatch?
 
 % % Both:
 %   residual2 = sqrt( mean( ( (abs(Hat_velocity) - abs(surf_vel_estimate2))/std_dev ).^2 ) ) + ...
 %            sqrt( mean( ( (abs(Hat_surface) - abs(S_P2(1,:)))/std_dev ).^2 ) );    
   residual2 = sqrt( mean( ( (abs(Hat_velocity(xpos2)) - abs(surf_vel_estimate2(xpos2)))/std_dev ).^2 ) ) + ...
              sqrt( mean( ( (abs(Hat_surface(xpos2)) - abs(S_P2(1,xpos2)))/std_dev ).^2 ) );   
          
  
  
  
 if (isreal(S_P2(time2,:)) == 0)
   disp('Imaginary in surface values')
   residual2 = NaN;   % Don't pick from this value if gives imaginary.
 end
  
 
   RMS_mismatch_matrix2(i_E,i_fs) = residual2;
   
   
   end % loop over fs values
   
  
  end  % loop over E values
    
    
   RMS_mismatch_matrix2;
   min_value2 = min(min(RMS_mismatch_matrix2));
  [minRow2, minCol2] = find(RMS_mismatch_matrix2 == min_value2);
 
   
%    if (length(index)>1)  % Multiple values of E give same mismatch
%        disp('Multiple values of fs give same mismatch value -- set fs=1e-11')
%        index = find(fs_vec == 1e-11);
%    end


   RMS_best               = RMS_mismatch_matrix2(minRow2, minCol2);
   best_fs                = fs_vec(minCol2);
   best_E                 = E_vec(minRow2);
   save_fs2(xpos2)          = best_fs;   
   save_E2(xpos2)           = best_E;
   save_mismatch2(xpos2, :) = RMS_mismatch_matrix2(minRow2, minCol2);

   
 % Keep best value:
   fs_P2(xpos2) = best_fs;                           
   E_P2(xpos2)  = best_E;
 
     
 
end  % loop over xpositions


      x_P2_min = x_P2;
      x_w2_min = x_w2;
      x_e2_min = x_e2;
      
      % Get around issue that interpolating to edge points can give negative  
      save_E2(1)   = save_E2(2);
      save_E2(end) = save_E2(end-1);
      E_P2_min     = save_E2;
      
    % recalculate values on the edges.
      [ E_w2_min, E_e2_min ] = get_edge_values_quadratic ...
                               ( E_P2, x_P2, x_w2, x_e2, dx_P2, dx_w2, dx_e2 );

     save_fs2(1)   = save_fs2(2);
     save_fs2(end) = save_fs2(end-1);
     fs_P2_min     = save_fs2;

    % recalculate values on the edges.
      [ fs_w2_min, fs_e2_min ] = get_edge_values_quadratic ...
                               ( fs_P2, x_P2, x_w2, x_e2, dx_P2, dx_w2, dx_e2 );
                       
   % Turns out that negatives can still come in. Just reset.                          
   fs_w2_min(1)   = fs_w2_min(2);
   fs_e2_min(end) = fs_e2_min(end-1);
   E_w2_min(1)    = E_w2_min(2);
   E_e2_min(end)  = E_e2_min(end-1);    
                           
   
    save min_Hat_E_and_fs.mat E_P2_min E_w2_min E_e2_min ...
                          fs_P2_min fs_w2_min fs_e2_min ...
                          x_P2_min x_w2_min x_e2_min

                           
  
                      
                     
  [ S_P2(time2,:), h_P2(time,:), ...  
    dS_dx_P_xt2(time2,:),  ...
    dS_dx_edges_xt2(time2,:) ] = calc_h_0( x_P2, x_w2, x_e2, dx_P2, B_P2, B_w2, B_e2, ...
                                         W_P2, W_w2, W_e2, b_dot_edges2(time2,:), ...
                                         E_w2, E_e2, fs_w2, fs_e2, ...
                                         Q_0_in2, S_at_GL2(time2), ...
                                         A_eff_edges_xt2(time2,1:end-1), ...
                                         A_eff_edges_xt2(time2,2:end), ...
                                         flux_add_w2, flux_add_e2, ...
                                         deformation_only2, deformation_plus_sliding2, sliding_only2);
              
% Find the flux
% ==============
 flux_kin_P_xt2(time2,:) = calc_flux_kin( x_P2, x_P2, dx_P2, W_P2, ...
                                        h_dot2(time2,:), b_dot_P2(time2,:), Q_0_in2 );
                                                                                             
[ flux_w2, flux_e2 ] = get_edge_values_quadratic ( flux_kin_P_xt2(time2,:), ...
                                                 x_P2, x_w2, x_e2, ...
                                                 dx_P2, dx_w2, dx_e2 );
   
 flux_edges_kin_xt2(time2,:) = [ Q_0_in2 flux_e2];   % Need to replace with Q_0_in, not extrapolated value.
                            
         
 % Velocity
 % ========
 surf_vel_estimate2 = (5/4) * abs(flux_edges_kin_xt2(1,:)) ./ ([W_w2(1) W_e2] .* [h_w2(1,1) h_e2(1,:)]);

    
 
 
 

figure(11)   % SURFACE and BED -- Hatherton
subplot(2,1,1), plot(x_P2/1000, B_P2,'r')
hold on
plot(x_P2/1000, S_modern2, 'linewidth', 2)
plot(x_P2/1000, S_P2(1,:),'c--')
legend('Bed', 'Measured surface', 'Calculated surface', 'location', 'northwest')
title('Surface and Bed topography: HATHERTON')
%xlabel('Distance along flowband (km)')
ylabel('Elevation (m)')
xlim([x_P2(1) x_P2(end)]/1000)

subplot(2,1,2)
plot(x_P2/1000, S_modern2-S_P2(1,:), 'linewidth', 2)
hold on
plot([x_P2(1) x_P2(end)]/1000, [0 0],'k--')
title('Modern surface - calculated surface')
xlabel('Distance along flowband (km)')
ylabel('Elevation difference (m)')
xlim([x_P2(1) x_P2(end)]/1000)



figure(31)   % SURFACE VELOCITY 
subplot(2,1,1), plot(x_edges2/1000, surf_vel_estimate2*(5/4), 'k--', 'linewidth', 2)
hold on
plot(Hat_measures_centerline_distance/1000, Hat_measures_flowspeed,'m', 'linewidth', 2)
title('Surface velocity: HATHERTON')
ylabel('Velocity (m/yr)')
legend('Calculated velocity',  'MEaSUREs velocity', 'location', 'best')

measures_on_xedges2 = interp1( Hat_measures_centerline_distance, Hat_measures_flowspeed, x_edges2, 'linear', 'extrap');

subplot(2,1,2), plot(x_edges2/1000, measures_on_xedges2 - average_vel_estimate2*(5/4), 'b', 'linewidth', 2)
hold on
plot([x_edges2(1) x_edges2(end)]/1000, [0 0],'k--')
title('Observed - calculated surface velocity')
xlabel('Distance along flowband (km)')
ylabel('Velocity difference (m/yr)')


                      
    end                        
 
 end
 
 
                                                  
      