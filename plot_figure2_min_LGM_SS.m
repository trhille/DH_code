
% ---------------
% Figure to plot finding LGM elevation at Darwin that best matches
% Hatherton data
% ---------------

% a) RMS mismatch vs. LGM elevation vector
% b) Darwin surfaces for LGM SS 
% c) Hatherton mismatch to data 
% d) Hatherton surfaces for LGM SS 

addpath('export_fig')

run_value = [1 2 4:7];


for jj = 1:length(run_value)
ii = run_value(jj);
 eval(['load LGM_steady_state_min' int2str(ii) '.mat'])
 
    
figure(2)
set(gcf, 'Units', 'centimeters','position', [35 20 24 20])
subplot('position', [0.15 0.61 0.2 0.34])
set(gca, 'fontsize', 16) 
% mismatch vs. LGM elevation values
plot(LGM_elev_vec, RMS_mismatch, 'color', [0.7 0.7 0.7], 'linewidth', 1)
hold on
set (gca, 'fontsize', 16)
xlabel({'LGM elevation ';'at Darwin Gl. (m)'}, 'fontsize', 16, 'fontweight', 'bold')
ylabel({'RMS mismatch'; 'along Hatherton Gl.'}, 'fontsize', 16, 'fontweight', 'bold')
if (ii == 1)
 plot(LGM_elev_vec, RMS_mismatch, 'k', 'linewidth', 2)
end    
xlim([300 1400])
set(gca,'xtick', [400 900 1400],'xticklabel',{'400','900','1400'})

subplot('position', [0.5 0.61 0.45 0.34])
set(gca, 'fontsize', 16) 
plot(x_P/1000, B_P,'r')
hold on
plot(x_P/1000, S_modern, 'linewidth', 2)
plot(x_P/1000, S_P(1,:),'c')
if (ii == 1)
h = legend('Bed', 'Modern surface', 'Calculated LGM', 'location', 'southeast');
set(h, 'fontsize', 12)
plot(x_P/1000, S_P(1,:),'k', 'linewidth', 1)
text(18, 2100, 'Darwin Glacier', 'fontangle', 'italic', 'fontsize', 16)
end
ylabel('Elevation (m)', 'fontweight', 'bold')
xlabel('Distance along flowband (km)', 'fontweight', 'bold')
xlim([x_P(1) x_P(end)]/1000)
ylim([-1000 2400])


   DAN_calc_elevation = interp1( x_P2, S_P2(time,:), DAN_position ); 
   MVfloor_calc_elevation = interp1( x_P2, S_P2(time,:), MVfloor_position ); 
   LWC14_calc_elevation = interp1( x_P2, S_P2(time,:), LWC14_position ); 
   
subplot('position', [0.15 0.12 0.2 0.34])
set(gca, 'fontsize', 16) 
plot(DAN_position/1000, DAN_mean_elevation - DAN_calc_elevation, 'g.', 'markersize', 12)
hold on
plot(MVfloor_position/1000, MVfloor_mean_elevation - MVfloor_calc_elevation, 'm.', 'markersize', 12)
plot(LWC14_position/1000, LWC14_mean_elevation - LWC14_calc_elevation, 'k.', 'markersize', 12)
if (ii == 1)
    plot(DAN_position/1000, DAN_mean_elevation - DAN_calc_elevation, 'o', 'markersize', 10,'markeredgecolor', 'k', 'markerfacecolor', 'g')
    plot(MVfloor_position/1000, MVfloor_mean_elevation - MVfloor_calc_elevation, 'o', 'markersize', 10,'markeredgecolor', 'k', 'markerfacecolor', 'm')
    plot(LWC14_position/1000, LWC14_mean_elevation - LWC14_calc_elevation, 'o', 'markersize', 10,'markeredgecolor', 'k', 'markerfacecolor', 'k')
end
xlabel({'Position along'; 'flowband (km)'}, 'fontweight', 'bold')
ylabel({'Elevation';'difference (m)'}, 'fontweight', 'bold')



subplot('position', [0.5 0.12 0.45 0.34])
set(gca, 'fontsize', 16) 
plot(x_P2/1000, B_P2,'r')
hold on
plot(x_P2/1000, S_modern2, 'linewidth', 2)
plot(x_P2/1000, S_P2(1,:),'c')
plot(DAN_position/1000, DAN_mean_elevation, 'g*')
plot(MVfloor_position/1000, MVfloor_mean_elevation, 'm*')
plot(LWC14_position/1000, LWC14_mean_elevation, 'k*')
 plot([DAN_position DAN_position]/1000, [700 1500], 'g')
 plot([MVfloor_position MVfloor_position]/1000, [700 1500], 'm')
 plot([LWC14_position LWC14_position]/1000, [700 1500], 'k')
if (ii == 1)
 h = legend('Bed', 'Modern surface', 'Calculated LGM', 'location', 'southeast');
set(h, 'fontsize', 12)
plot(x_P2/1000, S_P2(1,:),'k', 'linewidth', 1)
text(14, 2100, 'Hatherton Glacier', 'fontangle', 'italic', 'fontsize', 16)
 text(11, 1650, 'LW', 'fontsize', 14)
 text(28, 1650, 'MV', 'fontsize', 14, 'color', 'm')
 text(42, 1650, 'DAN', 'fontsize', 14, 'color', 'g')
end
xlabel('Distance along flowband (km)', 'fontweight', 'bold')
ylabel('Elevation (m)', 'fontweight', 'bold')
xlim([x_P2(1) x_P2(end)]/1000)
ylim([-500 2300])


end


export_fig figure2_LGM_SS.pdf -pdf -transparent


