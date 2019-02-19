function [ b_dot_nodes, ...
           b_dot_P, ...
           b_dot_edges ] = load_b_dot( x_P, x_edges, t_P, ...
                                       x_nodes, t_nodes )

% -------------------------------------------------------------------------                                   
                                   
global N_t_nodes N_x_nodes
  
global DIRECTORY_data
global precip_at_sl_LGM
global lapse_LGM 

addpath(DIRECTORY_data)

scale_Taylor_Dome = 0;

bdot1 = 0;
bdot2 = 0;
bdot3 = 0;
bdot4 = 0;
bdot5 = 1;
bdot6 = 1;
bdot7 = 0;
bdot8 = 0;

disp (' ')
disp ('Check bdot prescribed in load_b_dot.m')
disp(' ')


% Interpolate onto nodes and mesh:
% ================================
  b_dot_nodes = NaN * ones(N_t_nodes, N_x_nodes);
  
  
  
% Estimate from available observations
% ====================================

% From QGIS -- Accumulation_A vs. Accumulation_R
load DH_accum_width_velocity.mat
load DH_surf_bed.mat

%%

if bdot1 + bdot2 + bdot3 + bdot4 + bdot5 + bdot6 + bdot7 + bdot8 == 1

if (bdot1 == 1)
% % ------------------------------
% % OPTION 1: lapse rate estimate modern
% % ------------------------------

% if you want to calculate from a lapse rate similar to Bliss et al. (2011) for Taylor glacier 
precip_at_sl = -0.35;
lapse        = 0.35/1500;
Darwin_bdot_modern_lapse = precip_at_sl + lapse.*Darwin_modern_surface;
b_dot_use    = interp1(Darwin_centerline_distance, Darwin_bdot_modern_lapse, x_nodes);

elseif (bdot2 == 1)
% ------------------------------
% OPTION 2: Arthern estimate LGM
% ------------------------------
 b_dot_use = interp1(Darwin_accumulation_centerline_distance, Darwin_accumulation_A*0.6, x_nodes);  % NOT NEGATIVE here...


 elseif (bdot3 == 1)
% ------------------------------
% OPTION 3: RACMO estimate from Van de Berg et al. (2006)
% ------------------------------
 b_dot_use = interp1(Darwin_accumulation_centerline_distance, Darwin_accumulation_R, x_nodes);  % NOT NEGATIVE here...


 elseif (bdot4 == 1)
% ------------------------------
% OPTION 4: LGM estimate
% ------------------------------
  b_dot_use = interp1(Darwin_accumulation_centerline_distance, Darwin_accumulation_LGM, x_nodes);

elseif (bdot5 == 1)
% ------------------------------
% OPTION 5: Use RACMO2.1 5.5 km product
% ------------------------------
  load Darwin_RACMO2_1.mat
  
    b_dot_use = interp1(Darwin_distance_along_centerline + x_nodes(1), Darwin_SMB,...
      x_nodes, 'linear', 'extrap');

elseif (bdot6 == 1)
% ------------------------------
% OPTION 6: Use RACMO2.3 27 km product
% ------------------------------
  load Darwin_RACMO2_3.mat
  
    b_dot_use = interp1(Darwin_distance_along_centerline + x_nodes(1), Darwin_SMB,...
      x_nodes, 'linear', 'extrap');

end

  for ii = 1:N_t_nodes
    b_dot_nodes(ii,:) = b_dot_use;
  end
  
  %% Give n_spinup timesteps to equilibrate to new SMB if bdot5 ~= 1
    load Darwin_RACMO2_1.mat
    n_spinup = 10; %numer of nodes for spinup to account for different SMB than that used for min_E and min_fs (bdot5)
    b_dot_spinup = interp1(Darwin_distance_along_centerline + x_nodes(1), Darwin_SMB,...
      x_nodes, 'linear', 'extrap');

  spinup_weight = linspace(0,1,n_spinup);
  
  for mm = 1:n_spinup
  b_dot_nodes(mm,:) = spinup_weight(mm).*b_dot_use + ...
      (1-spinup_weight(mm)).*b_dot_spinup;
  end

elseif bdot1 + bdot2 + bdot3 + bdot4 + bdot5 + bdot6 + bdot7 + bdot8  > 1

if bdot1 == 1 && bdot2 == 1
    weight = linspace(1,0, N_t_nodes);
    
    
    precip_at_sl = -0.35;
    lapse = 0.35/1500;
    Darwin_bdot_modern_lapse = precip_at_sl + lapse.*Darwin_modern_surface;
    b_dot_modern = interp1(Darwin_centerline_distance, Darwin_bdot_modern_lapse, x_nodes);
    b_dot_LGM = interp1(Darwin_accumulation_centerline_distance, Darwin_accumulation_A*0.6, x_nodes);

    for ii = 1:N_t_nodes
        b_dot_nodes(ii,:) = b_dot_LGM.*weight(ii) + b_dot_modern.*(1-weight(ii));
    end
    
    
elseif bdot5 == 1 && bdot6 == 1
    load Darwin_RACMO2_3.mat
    Darwin_distance_along_centerline2_3 = Darwin_distance_along_centerline;
    Darwin_SMB2_3 = Darwin_SMB;
    b_dot_LGM = 0.6.*interp1(Darwin_distance_along_centerline2_3 + x_nodes(1), Darwin_SMB2_3,...
        x_nodes, 'linear', 'extrap');   
    
    load Darwin_RACMO2_1.mat
    Darwin_distance_along_centerline2_1 = Darwin_distance_along_centerline;
    Darwin_SMB2_1 = Darwin_SMB; clear Darwin_distance_along_centerline Darwin_SMB
    b_dot_modern = interp1(Darwin_distance_along_centerline2_1 + x_nodes(1), Darwin_SMB2_1,...
        x_nodes, 'linear', 'extrap');

    weight = [linspace(0,1, ceil(N_t_nodes ./4)), ones(1, ceil(N_t_nodes ./2)), linspace(1,0, ceil(N_t_nodes ./4))];
    
    for ii = 1:N_t_nodes
        b_dot_nodes(ii,:) = b_dot_LGM.*weight(ii) + b_dot_modern.*(1-weight(ii));
    end

elseif bdot5 == 1 && bdot7 == 1
load Darwin_RACMO2_3.mat
Darwin_distance_along_centerline2_3 = Darwin_distance_along_centerline;
Darwin_SMB2_3 = Darwin_SMB;
b_dot_LGM = interp1(Darwin_distance_along_centerline2_3 + x_nodes(1), Darwin_SMB2_3,...
    x_nodes, 'linear', 'extrap');   

load Darwin_RACMO2_1.mat
Darwin_distance_along_centerline2_1 = Darwin_distance_along_centerline;
Darwin_SMB2_1 = Darwin_SMB; clear Darwin_distance_along_centerline Darwin_SMB
b_dot_modern = interp1(Darwin_distance_along_centerline2_1 + x_nodes(1), Darwin_SMB2_1,...
    x_nodes, 'linear', 'extrap');

weight = [linspace(0,1, ceil(N_t_nodes ./4)), ones(1, ceil(N_t_nodes ./2)), linspace(1,0, ceil(N_t_nodes ./4))];

for ii = 1:N_t_nodes
    b_dot_nodes(ii,:) = b_dot_LGM.*weight(ii) + b_dot_modern.*(1-weight(ii));
end

elseif bdot5 == 1 && bdot8 == 1
load Darwin_RACMO2_3.mat
Darwin_distance_along_centerline2_3 = Darwin_distance_along_centerline;
Darwin_SMB2_3 = Darwin_SMB;
b_dot_LGM = 2*interp1(Darwin_distance_along_centerline2_3 + x_nodes(1), Darwin_SMB2_3,...
    x_nodes, 'linear', 'extrap');   

load Darwin_RACMO2_1.mat
Darwin_distance_along_centerline2_1 = Darwin_distance_along_centerline;
Darwin_SMB2_1 = Darwin_SMB; clear Darwin_distance_along_centerline Darwin_SMB
b_dot_modern = interp1(Darwin_distance_along_centerline2_1 + x_nodes(1), Darwin_SMB2_1,...
    x_nodes, 'linear', 'extrap');

weight = [linspace(0,1, ceil(N_t_nodes ./4)), ones(1, ceil(N_t_nodes ./2)), linspace(1,0, ceil(N_t_nodes ./4))];

for ii = 1:N_t_nodes
    b_dot_nodes(ii,:) = b_dot_LGM.*weight(ii) + b_dot_modern.*(1-weight(ii));
end    
end

end 
%% This scales a given modern accumulation rate to the Taylor Dome record 
%from Monnin et al. (2004) in terms of % of modern accumulation rate. This
%can be used with any combination of bdotx above, but only really makes
%sense to use bdot5 or bdot6

if scale_Taylor_Dome == 1
    TD_record = readtable('DH_DATA/Taylor_Dome_accum_rate.csv');
    
     % use nearest interpolation to pin the record at the edges so it
     % doesn't grow enormous towards the modern
    TD_accum_rate_interp = interp1(TD_record.yr_BP, ...
        TD_record.accum_rate, t_nodes, 'nearest', 'extrap');
    
    %normalize to modern
    TD_accum_frac_modern = smooth(TD_accum_rate_interp./TD_accum_rate_interp(end));
    
    n_t_spinup = 15; %number of timesteps used for spinup to new accum rate
    TD_accum_frac_modern(1:n_t_spinup) = linspace(1, ...
        TD_accum_frac_modern(n_t_spinup), n_t_spinup);
     
    for ii = 1:N_t_nodes
    b_dot_nodes(ii,:) = b_dot_nodes(ii,:) .* TD_accum_frac_modern(ii);
    end
end  
    


  
 %% 11/25/18 Experimenting with varying accumulation in time
% load('TEST_DH/output/LGM_b_dot/precip_0_2lapse_-6_6667e-05.mat', 'b_dot_P', 'b_dot_edges');
% b_dot_P_LGM = b_dot_P(end,:); b_dot_edges_LGM = b_dot_edges; clear b_dot_P b_dot_edges;
% 
% b_dot_P_modern = interp1(x_nodes, b_dot_use, x_P);
% 
% weight = linspace(1,0, length(t_P));
% 
% for ii = 1:length(t_P)
%     b_dot_P(ii,:) = b_dot_P_LGM.*weight(ii) + b_dot_P_modern.*(1-weight(ii));
% end
% 
% % Constant in time
% % ------------------
%  b_dot_use = interp1([x_P(1) x_P(end)], [0.3 0.3], x_P);  % Average accumulation 28 cm/yr
%                                                           % Denton et al. (1989)                                                                                                          
%  for ii = 1:N_t_nodes
%    b_dot_nodes(ii,:) = interp1(x_P, b_dot_use, ...
%                                x_nodes, 'linear', 'extrap');
%  end
% 
%        
% % Linearly decreasing in time
% % ----------------------------
% b_dot_use      = interp1([x_P(1) x_P(end)], [0.1 0.1], x_P);
% time_variation = interp1([t_P(1) t_P(end)], [6 3], t_nodes); 
% 
% for ii = 1:N_t_nodes
%     b_dot_nodes(ii,:) = time_variation(ii) * interp1(x_P, b_dot_use, x_nodes, 'linear', 'extrap');
% end




% interpolate onto mesh -- from nodes!
% ----------------------                                        
   b_dot_P = interp2( repmat( x_nodes, N_t_nodes, 1 ), ...
                      repmat( t_nodes, 1, N_x_nodes ), ...
                      b_dot_nodes, x_P, t_P );

 
% interpolate onto edges of mesh
% ------------------------------ 
   b_dot_edges = interp2( repmat( x_nodes, N_t_nodes, 1 ), ...
                           repmat( t_nodes, 1, N_x_nodes ), ...
                           b_dot_nodes, x_edges, t_P );  
                       
%      b_dot_edges = interp2( repmat( x_P, length(t_P), 1 ), ...
%                            repmat( t_P, 1, length( x_P) ), ...
%                            b_dot_P, x_edges, t_P );   
%                        b_dot_edges(:,1) = b_dot_edges(:,2);
end                                              