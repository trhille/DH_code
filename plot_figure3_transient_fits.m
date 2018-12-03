
% -------------------------------------------------------------------------
% Results from transient runs -- prescribe LGM elevation history at Darwin;
% compare to data from Hatherton
% -------------------------------------------------------------------------

addpath('export_fig')

load LGM_steady_state_min2.mat   % to load mean elevation values of data.


  DIRECTORY_data = 'DH_DATA';
  addpath( DIRECTORY_data )


% NOTE: Need to change file names in the plot command below for smooth or step!


run_value = [1:7];  

for kk = 1:length(run_value)
   
    jj = run_value(kk);


 eval(['load smooth_9ka_min' int2str(jj) '.mat'])
% eval(['load stepwise_9ka_min' int2str(jj) '.mat'])
%   eval(['load pollard_min' int2str(jj) '.mat'])


 index_val = jj;
 
figure
set(gcf, 'Units', 'centimeters','position', [35 20 24 30])
subplot('position', [0.15 0.875 0.8 0.1])
    plot(-t_P/1000, S_at_GL - S_at_GL(end),'k', 'linewidth', 2)
    set(gca, 'fontsize', 16)
    xlabel('Time (kyr)', 'fontweight', 'bold')
    ylabel({'Height above'; 'modern (m)'}, 'fontweight', 'bold')
    %text(-11.6, 400, 'Prescribed elevation near ', 'fontangle', 'italic', 'fontsize', 14)
    %text(-11.6, 200, 'Darwin Glacier grounding line', 'fontangle', 'italic', 'fontsize', 14)
    xlim([0 12])
    ylim([-50 700]) 
    
subplot('position', [0.15 0.61 0.8 0.2])
    plot(x_P/1000, S_P(81:20:end,:),'b')
    hold on
    plot(x_P/1000, B_P,'r')
    plot([x_P(index_xpos_Hatherton_P) x_P(index_xpos_Hatherton_P)]/1000,[B_P(index_xpos_Hatherton_P) 2200], 'c')
    plot(x_P/1000, S_modern(end,:),'k--','linewidth', 2)
    set(gca, 'fontsize', 16)
    %xlabel('Distance along flowband (km)', 'fontweight', 'bold')
    ylabel({'Elevation'; '(m a.s.l.)'}, 'fontweight', 'bold')
    xlim([17 150])
    ylim([-750 2200])
    text(120, -300, 'Darwin Glacier', 'fontangle', 'italic', 'fontsize', 18)
    text(36, 1850, 'Hatherton Glacier', 'fontangle', 'italic', 'fontsize', 14)
    text(36, 1500, 'input', 'fontangle', 'italic', 'fontsize', 14)

    
 load '/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton/Modeling/Koutnik model /DH_code/DH_DATA/all_values.mat'   
 load '/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton/Modeling/Koutnik model /DH_code/DH_DATA/LGM_values_use.mat'
 index_DAN     = find(x_P2 <= DAN_position, 1, 'last');
 index_MVfloor = find(x_P2 <= MVfloor_position, 1, 'last')-1;
 index_LWC14   = find(x_P2 <= LWC14_position, 1, 'last'); 
    
 
subplot('position', [0.15 0.36 0.8 0.2])
    plot(x_P2/1000, S_P2(81:20:end,:),'b')
    hold on
    plot(x_P2/1000, B_P2,'r')
    plot(x_P2/1000, S_modern2, 'k--', 'linewidth', 2)
    set(gca, 'fontsize', 16)
    xlabel('Distance along flowband (km)', 'fontweight', 'bold')
    ylabel({'Elevation'; '(m a.s.l.)'}, 'fontweight', 'bold')
    ylim([-200 2100])
    text(65, 100, 'Hatherton Glacier', 'fontangle', 'italic', 'fontsize', 18)
    xlim([10 85])
%     plot([DAN_position DAN_position]/1000, [700 1500], 'g')
%     plot([MVfloor_position MVfloor_position]/1000, [700 1500], 'm')
%     plot([LWC14_position LWC14_position]/1000, [700 1500], 'k')
    plot(DAN_position/1000, DAN_mean_elevation, 'g*')
    plot(MVfloor_position/1000, MVfloor_mean_elevation, 'm*')
    plot(LWC14_position/1000, LWC14_mean_elevation, 'k*')
    text(11, 1700, 'LW', 'fontsize', 14)
    text(28, 1700, 'MV', 'fontsize', 14, 'color', 'm')
    text(43, 1700, 'DAN', 'fontsize', 14, 'color', 'g')
   
    
 for iii = 1:length(run_value)
 
     ii = run_value(iii);
     
 eval(['load smooth_9ka_min' int2str(ii) '.mat'])
 % eval(['load stepwise_9ka_min' int2str(ii) '.mat'])     
%    eval(['load pollard_min' int2str(ii) '.mat'])

  subplot('position', [0.15 0.075 0.2 0.2])
    plot(-t_P2/1000, S_P2(:,index_LWC14), 'color', [0.7 0.7 0.7], 'linewidth', 2)
    hold on
    plot(LWC14_ages/1000, LWC14_elevations, 'k.', 'markersize', 12)
    xlim([0 12])
    set(gca, 'fontsize', 16)
    ylabel({'Elevation'; '(m a.s.l.)'}, 'fontweight', 'bold')
    
 subplot('position', [0.45 0.075 0.2 0.2])
    plot(-t_P2/1000, S_P2(:,index_MVfloor), 'color', [0.7 0.7 0.7], 'linewidth', 2)
    hold on
    plot(MVfloor_Be10_ages/1000, MVfloor_elevations, 'm.', 'markersize', 12)
    xlim([0 12])
    set(gca, 'fontsize', 16)
    xlabel('Time (kyr)', 'fontweight', 'bold')
         
 subplot('position', [0.75 0.075 0.2 0.2])
    plot(-t_P2/1000, S_P2(:,index_DAN), 'color', [0.7 0.7 0.7], 'linewidth', 2)
    hold on
    plot(DAN_Be10_ages/1000, DAN_elevations, 'g.', 'markersize', 12)
    xlim([0 12])
    set(gca, 'fontsize', 16)
    
end

 eval(['load smooth_9ka_min' int2str(index_val) '.mat'])
%   eval(['load stepwise_9ka_min' int2str(index_val) '.mat'])
%   eval(['load pollard_min' int2str(index_val) '.mat'])


 subplot('position', [0.15 0.075 0.2 0.2])
    plot(-t_P2/1000, S_P2(:,index_LWC14), 'c', 'linewidth', 2)
    hold on
    plot(LWC14_ages/1000, LWC14_elevations, 'k.', 'markersize', 12)
    xlim([1 12])
    set(gca, 'fontsize', 16)
    ylabel({'Elevation'; '(m a.s.l.)'}, 'fontweight', 'bold')
    
 subplot('position', [0.45 0.075 0.2 0.2])
    plot(-t_P2/1000, S_P2(:,index_MVfloor), 'c', 'linewidth', 2)
    hold on
    plot(MVfloor_Be10_ages/1000, MVfloor_elevations, 'm.', 'markersize', 12)
    xlim([0 12])
    set(gca, 'fontsize', 16)
    xlabel('Time (kyr BP)', 'fontweight', 'bold')
         
 subplot('position', [0.75 0.075 0.2 0.2])
    plot(-t_P2/1000, S_P2(:,index_DAN), 'c', 'linewidth', 2)
    hold on
    plot(DAN_Be10_ages/1000, DAN_elevations, 'g.', 'markersize', 12)
    xlim([0 12])
    set(gca, 'fontsize', 16)
    

    
% eval(['export_fig figure3_LGM_transient_smooth9ka_min' int2str(index_val) '.pdf -pdf -transparent'])
% eval(['export_fig figure3_LGM_transient_step9ka_min' int2str(index_val) '.pdf -pdf -transparent'])
%  eval(['export_fig figure3_LGM_transient_pollard_min' int2str(index_val) '.pdf -pdf -transparent'])

% close all

end


    
