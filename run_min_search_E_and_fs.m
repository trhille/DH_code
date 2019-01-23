% run_min_search_E


% % User needs to CHANGE iterations = 1 for Darwin, or = 2 for Haterhton!

iterations = 1   % Minimize for both Darwin (1) and Hatherton (2)
% iterations = 2    


E_vec         = [ 0.1 0.5 1 5 10 ];
%E_vec = [1]; % 3 10];
% fs_vec        = [ 0.1e-11 0.2e-11 0.3e-11 0.4e-11 0.5e-11 0.6e-11 0.7e-11 0.8e-11 0.9e-11 1e-11];
fs_vec = [ 1e-13 1e-12 1e-11];


N_E           = length( E_vec );
N_fs          = length( fs_vec );

% Still find best parameter set for Darwin, then for Hatherton
% But now search over combinations of E and sliding, include E=1 and fs=0
% Evaluate on edge values so don't get any negatives interpolating from
% centerpoints to edges

    
   
if (iterations == 1)
        disp('Running for Darwin...')
        disp(' ')
        length_run          = length(x_edges);
        save_mismatch       = NaN * ones(length_run, N_E*N_fs);
        E_edges_orig        = [E_w(1) E_e];
        fs_edges_orig       = [fs_w(1) fs_e];
        Darwin_velocity     = interp1(Darwin_measures_centerline_distance, ...
                                      Darwin_measures_flowspeed, [x_w(1) x_e]);
        Darwin_surface      = S_modern;
       [ Darwin_S_w ...
         Darwin_S_e]        = get_edge_values_quadratic ( S_modern, ...
                                                          x_P, x_w, x_e, ...
                                                          dx_P, dx_w, dx_e );
        Darwin_surface_edges = [ Darwin_S_w(1) Darwin_S_e ];
        
        RMS_mismatch_matrix = NaN * ones(N_E, N_fs);
        
        E_edges_temp        = [E_w(1) E_e];
        fs_edges_temp       = [fs_w(1) fs_e];
        

 for xpos = 1:length_run   % loop over all x positions -- upstream
    
disp(['Evaluating E and fs for x-position ',int2str(xpos), ' of total=',int2str(length_run)]);

RMS_mismatch_matrix = NaN * ones( N_E, N_fs );


%  Loop over E_x
 for i_E = 1:N_E
     
      E_use              = E_vec(i_E);
      E_edges_temp(xpos) = E_use;
      
      E_w = E_edges_temp(1:end-1);
      E_e = E_edges_temp(2:end);
      E_P = interp1(x_edges, E_edges_temp, x_P);
  
                                              
   for i_fs = 1:N_fs                       
           
      fs_use               = fs_vec(i_fs);
      fs_edges_temp(xpos)  = fs_use;
    
      fs_w = fs_edges_temp(1:end-1);
      fs_e = fs_edges_temp(2:end);
      fs_P = interp1(x_edges, fs_edges_temp, x_P);
     
     
                           
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
                   
 
 
 
%  if (flux_edges_kin_xt(time, end) ~= 0)
%      
%    Q_0_in = Q_0_in - (flux_edges_kin_xt(time,end));
%    
%   [ S_P(time,:), h_P(time,:), ...  
%     dS_dx_P_xt(time,:),  ...
%     dS_dx_edges_xt(time,:) ] = calc_h_0( x_P, x_w, x_e, dx_P, B_P, B_w, B_e, ...
%                                          W_P, W_w, W_e, b_dot_edges(time,:), ...
%                                          E_w, E_e, fs_w, fs_e, ...
%                                          Q_0_in, S_at_GL(time), ...
%                                          A_eff_edges_xt(time,1:end-1), ...
%                                          A_eff_edges_xt(time,2:end), ...
%                                          flux_add_w, flux_add_e, ...
%                                          deformation_only, deformation_plus_sliding, sliding_only);
%                                                               
%   [ h_w(time,:), ...
%     h_e(time,:) ] = get_edge_values_quadratic( h_P(time,:), x_P, x_w, x_e, ...
%                                                dx_P, dx_w, dx_e );
%  
%                                            
%  end                                          
 
 
 flux_kin_P_xt(time,:) = calc_flux_kin( x_P, x_P, dx_P, W_P, ...
                                h_dot(time,:), (b_dot_P(time,:)), Q_0_in );
                                                
 
 
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
   surf_vel_estimate = (5/4) * abs(flux_edges_kin_xt(1,:)) ./ ...
                       ([W_w(1) W_e] .* [h_w(1,1) h_e(1,:)]);  % USING KIN FLUX HERE!
 
                 
   [ S_w, S_e ] = get_edge_values_quadratic ( S_P(time,:), ...
                                                 x_P, x_w, x_e, ...
                                                 dx_P, dx_w, dx_e );
               
   S_edges = [ S_w(1,1) S_e(1,:) ];
       
 % Residual
 % ========
    std_dev = 1; 

%  % Surface profile:
%   residual = sqrt( mean( ( (abs(Darwin_surface) - abs(S_P(1,:)))/std_dev ).^2 ) );           
%   residual = sqrt( mean( ( (abs(Darwin_surface(xpos)) - abs(S_P(1,xpos)))/std_dev ).^2 ) );  

%  % Surface velocity:
%   residual = sqrt( mean( ( (abs(Darwin_velocity) - abs(surf_vel_estimate))/std_dev ).^2 ) );  
%   residual = abs(Darwin_velocity(xpos) - surf_vel_estimate(xpos));  % local mismatch?
 
% % Both:
%    residual = sqrt( mean( ( (abs(Darwin_velocity) - abs(surf_vel_estimate))/std_dev ).^2 ) ) + ...
%               sqrt( mean( ( (abs(Darwin_surface_edges) - abs(S_edges(1,:)))/std_dev ).^2 ) );   

  residual = sqrt( mean( ( (abs(Darwin_velocity(xpos)) - abs(surf_vel_estimate(xpos)))/std_dev ).^2 ) ) + ...
             sqrt( mean( ( (abs(Darwin_surface_edges(xpos)) - abs(S_edges(xpos)))/std_dev ).^2 ) );   
          
  
if (isreal(S_P(1,:)) == 0)
    residual = NaN;   % Don't pick from this value if gives imaginary.
end
  
   RMS_mismatch_matrix(i_E,i_fs) = residual; % this resets for each xpos, 
                                             % used only to find best value at xpos
   

   
   end % loop over fs values
   
  
  end  % loop over E values
    
    
   min_value = min(min(RMS_mismatch_matrix));
  [minRow, minCol] = find(RMS_mismatch_matrix == min_value);
   
  
    if (length(minRow)>1)  % Multiple values of E give same mismatch
        disp('Multiple values of E give same mismatch value -- line 198 run_min_search_E_and_fs')
        best_fs = fs_edges_temp(xpos-1); % set to be neighboring value
        
    elseif (length(minCol)>1) % Multiple values of fs give same mismatch
        disp('Multiple values of fs give same mismatch value -- line 198 run_min_search_E_and_fs')
        best_E = E_edges_temp(xpos-1);  % set to be neighboring value
          
    elseif (length(minCol)==1) && (length(minRow)==1 )
       RMS_best               = RMS_mismatch_matrix(minRow, minCol);
       best_fs                = fs_vec(minCol);
       best_E                 = E_vec(minRow);
       save_mismatch(xpos, :) = RMS_mismatch_matrix(minRow, minCol); 
    end % ifstatement on whether there is one minimum

    
 

   E_edges_temp(xpos)  = best_E;
   fs_edges_temp(xpos) = best_fs;
   

end  % loop over xpositions



%   E_w = E_edges_temp(1:end-1);
%   E_e = E_edges_temp(2:end);
%   E_P_min = interp1(x_edges, E_edges_temp, x_P);
%   E_w_min = interp1(x_P, E_P_min, x_edges(1:end-1), 'linear', 'extrap');
%   E_e_min = interp1(x_P, E_P_min, x_edges(2:end), 'linear', 'extrap');
%   
%   fs_w = fs_edges_temp(1:end-1);
%   fs_e = fs_edges_temp(2:end);
%   fs_P_min = interp1(x_edges, fs_edges_temp, x_P);
%   fs_w_min = interp1(x_P, fs_P_min, x_edges(1:end-1),'linear', 'extrap');
%   fs_e_min = interp1(x_P, fs_P_min, x_edges(2:end),'linear', 'extrap');

  E_w_min = E_edges_temp(1:end-1);
  E_e_min = E_edges_temp(2:end);
  E_P_min = NaN;
  
  fs_w_min = fs_edges_temp(1:end-1);
  fs_e_min = fs_edges_temp(2:end);
  fs_P_min = NaN;
  
  x_P_min = x_P;
  x_w_min = x_w;
  x_e_min = x_e;


save min_Darwin_E_and_fs.mat E_P_min E_w_min E_e_min ...
                             fs_P_min fs_w_min fs_e_min ...
                             x_P_min x_w_min x_e_min  
end  % loop on Darwin (1) or Hatherton (2)

 


 % ------------------------------------------------------------------------
 % ------------------------------------------------------------------------
 % ------------------------------------------------------------------------
 % ------------------------------------------------------------------------                          
 
 
if (iterations == 2)
        
        disp('Running for Hatherton...')
        disp(' ')
        length_run2          = length(x_edges2);
        save_mismatch2       = NaN * ones(length_run2, N_E);
        E_edges_orig2        = [E_w2(1) E_e2];
        fs_edges_orig2       = [fs_w2(1) fs_e2]; 
        Hat_velocity         = interp1( Hat_measures_centerline_distance, ...
                                        Hat_measures_flowspeed, x_edges2, ...
                                        'linear', 'extrap');
       [ Hat_S_w Hat_S_e]    = get_edge_values_quadratic ( S_modern2, ...
                                                           x_P2, x_w2, x_e2, ...
                                                           dx_P2, dx_w2, dx_e2 );
        Hat_surface_edges = [ Hat_S_w(1) Hat_S_e ];
       
        
        RMS_mismatch_matrix2 = NaN * ones(N_E, N_fs);
        
    
        E_edges_temp2        = [E_w2(1) E_e2];
        fs_edges_temp2       = [fs_w2(1) fs_e2];
        
        
 for xpos2 = 1:length_run2   % loop over all x edges positions -- upstream
    
disp(['Evaluating E and fs for x-position ',int2str(xpos2), ' of total=',int2str(length_run2)]);

RMS_mismatch_matrix2 = NaN * ones( N_E, N_fs );


%  Loop over E_x
 for i_E = 1:N_E
     
      E_use2               = E_vec(i_E);
      E_edges_temp2(xpos2) = E_use2;
      
      E_w2 = E_edges_temp2(1:end-1);
      E_e2 = E_edges_temp2(2:end);
      E_P2 = interp1(x_edges2, E_edges_temp2, x_P2);
  
                            
                           
   for i_fs = 1:N_fs                       
           
      fs_use2               = fs_vec(i_fs);
      fs_edges_temp2(xpos2) = fs_use2;
    
      fs_w2 = fs_edges_temp2(1:end-1);
      fs_e2 = fs_edges_temp2(2:end);
      fs_P2 = interp1(x_edges2, fs_edges_temp2, x_P2);

                             
                         
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
                            
 
 
%  if (flux_edges_kin_xt2(time2,end) ~= 0)
%   
%    Q_0_in2 = Q_0_in2 - (flux_edges_kin_xt2(time2, end));
%      
%    [ S_P2(time2,:), h_P2(time,:), ...  
%     dS_dx_P_xt2(time2,:),  ...
%     dS_dx_edges_xt2(time2,:) ] = calc_h_0( x_P2, x_w2, x_e2, dx_P2, B_P2, B_w2, B_e2, ...
%                                          W_P2, W_w2, W_e2, b_dot_edges2(time2,:), ...
%                                          E_w2, E_e2, fs_w2, fs_e2, ...
%                                          Q_0_in2, S_at_GL2(time2), ...
%                                          A_eff_edges_xt2(time2,1:end-1), ...
%                                          A_eff_edges_xt2(time2,2:end), ...
%                                          flux_add_w2, flux_add_e2, ...
%                                          deformation_only2, deformation_plus_sliding2, sliding_only2);
%                                      
%   [ h_w2(time2,:), ...
%     h_e2(time2,:) ] = get_edge_values_quadratic( h_P2(time2,:), x_P2, x_w2, x_e2, ...
%                                                  dx_P2, dx_w2, dx_e2 );
%                                                 
%  end
                                                                      
  
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
                                                                       
 flux_edges_dyn_xt2(1,1)    = Q_0_in2; 
                                        
                   
 [ S_w2, S_e2 ] = get_edge_values_quadratic ( S_P2(time,:), ...
                                              x_P2, x_w2, x_e2, ...
                                              dx_P2, dx_w2, dx_e2 );         
   S_edges2 = [ S_w2(1,1) S_e2(1,:) ];
 
 
 % Velocity
 % ========
 surf_vel_estimate2 = (5/4) * abs(flux_edges_kin_xt2(1,:)) ./ ([W_w2(1) W_e2] .* [h_w2(1,1) h_e2(1,:)]);
 average_vel_estimate2 = abs(flux_edges_kin_xt2(1,:)) ./ ([W_w2(1) W_e2] .* [h_w2(1,1) h_e2(1,:)]);  % USING KIN FLUX HERE!
 
 
 % Residual
 % ========
   std_dev = 1; 

% % Surface profile:
%   residual2 = sqrt( mean( ( (abs(Hat_surface) - abs(S_P2(1,:)))/std_dev ).^2 ) );           

% % Surface velocity:
%  residual2 = sqrt( mean( ( (abs(Hat_velocity) - abs(surf_vel_estimate2))/std_dev ).^2 ) );  
%  residual2 = abs(Hat_velocity(xpos2) - surf_vel_estimate2(xpos2));  % local mismatch?
 
 % % Both:
%  residual2 = sqrt( mean( ( (abs(Hat_velocity) - abs(surf_vel_estimate2))/std_dev ).^2 ) ) + ...
%               sqrt( mean( ( (abs(Hat_surface_edges) - abs(S_edges2(1,:)))/std_dev ).^2 ) );   
 
residual2 = sqrt( mean( ( (abs(Hat_velocity(xpos2)) - abs(surf_vel_estimate2(xpos2)))/std_dev ).^2 ) ) + ...
             sqrt( mean( ( (abs(Hat_surface_edges(xpos2)) - abs(S_edges2(xpos2)))/std_dev ).^2 ) );   
          
  
  
 if (isreal(S_P2(time2,:)) == 0)
   disp('Imaginary in surface values')
   residual2 = NaN;   % Don't pick from this value if gives imaginary.
 end
  
 
   RMS_mismatch_matrix2(i_E,i_fs) = residual2;
   
   
   
   end % loop over fs values
   
  
  end  % loop over E values
    
  
    
   min_value2 = min(min(RMS_mismatch_matrix2));
  [minRow2, minCol2] = find(RMS_mismatch_matrix2 == min_value2);
   
  
    if (length(minRow2)>1)  % Multiple values of E give same mismatch
        disp('Multiple values of E give same mismatch value -- line 198 run_min_search_E_and_fs')
        best_fs = fs_edges_temp2(xpos2-1); % set to be neighboring value
        
    elseif (length(minCol2)>1) % Multiple values of fs give same mismatch
        disp('Multiple values of fs give same mismatch value -- line 198 run_min_search_E_and_fs')
        best_E = E_edges_temp2(xpos2-1);  % set to be neighboring value
          
    elseif (length(minCol2)==1) && (length(minRow2)==1 )
       RMS_best               = RMS_mismatch_matrix2(minRow2, minCol2);
       best_fs                = fs_vec(minCol2);
       best_E                 = E_vec(minRow2);
       save_mismatch2(xpos2, :) = RMS_mismatch_matrix2(minRow2, minCol2); 
    end % ifstatement on whether there is one minimum

 
 
   E_edges_temp2(xpos2)  = best_E;
   fs_edges_temp2(xpos2) = best_fs;
  
      
 
end  % loop over xpositions


%   E_w2 = E_edges_temp2(1:end-1);
%   E_e2 = E_edges_temp2(2:end);
%   E_P2_min = interp1(x_edges2, E_edges_temp2, x_P2);
%   E_w2_min = interp1(x_P2, E_P2_min, x_edges2(1:end-1), 'linear', 'extrap');
%   E_e2_min = interp1(x_P2, E_P2_min, x_edges2(2:end), 'linear', 'extrap');
%   
%   fs_w2 = fs_edges_temp2(1:end-1);
%   fs_e2 = fs_edges_temp2(2:end);
%   fs_P2_min = interp1(x_edges2, fs_edges_temp2, x_P2);
%   fs_w2_min = interp1(x_P2, fs_P2_min, x_edges2(1:end-1),'linear', 'extrap');
%   fs_e2_min = interp1(x_P2, fs_P2_min, x_edges2(2:end),'linear', 'extrap');
  
  E_w2_min = E_edges_temp2(1:end-1);
  E_e2_min = E_edges_temp2(2:end);
  E_P2_min = NaN;
  
  fs_w2_min = fs_edges_temp2(1:end-1);
  fs_e2_min = fs_edges_temp2(2:end);
  fs_P2_min = NaN;

  x_P2_min = x_P2;
  x_w2_min = x_w2;
  x_e2_min = x_e2;

save min_Hat_E_and_fs.mat E_P2_min E_w2_min E_e2_min ...
                          fs_P2_min fs_w2_min fs_e2_min ...
                          x_P2_min x_w2_min x_e2_min

end  % loop over iterations 1 or 2
           



