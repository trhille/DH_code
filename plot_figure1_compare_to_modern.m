addpath('export_fig')


% ---------------
% Figure to show how these sets of BC match the modern data
% a) Darwin surface
% b) Darwin velocity
% c) Hatherton surface
% d) Hatherton velocity
% ---------------


run_value = [1 2 4:7];  

%run_value = 1

for ii = 1:length(run_value)
ii_use = run_value(ii);
%eval(['load run_min_search' int2str(ii_use) 'a.mat'])    
eval(['load min' int2str(ii_use) '.mat'])    

% Darwin : (modern surface - calculated surface)
figure(1)
set(gcf, 'Units', 'centimeters','position', [35 20 22 20])
subplot('position', [0.2 0.81 0.8 0.15])
set(gca, 'fontsize', 16)
plot(x_P/1000, (S_modern-S_P(1,:))./S_modern, 'color', [0.7 0.7 0.7], 'linewidth', 1)
hold on
plot([x_P(1) x_P(end)]/1000, [0 0],'k--')
if (ii == 1)
 plot(x_P/1000, (S_modern-S_P(1,:))./S_modern, 'k', 'linewidth', 2)
end    
ylabel({'Fractional';' elevation ';'difference ';'[obs-model]'}, 'fontweight', 'bold')
xlim([x_P(1) x_P(end)]/1000)
%ylim([-300 350])
ylim([-0.05 0.2])
if (ii == 6)
text(110, -100,'Darwin Glacier', 'fontangle', 'italic', 'fontsize', 16)
end
set(gca, 'xticklabel', [])


% Darwin : (observed - calculated surface velocity)
subplot('position', [0.2 0.6 0.8 0.15])
set(gca, 'fontsize', 16)
measures_on_xedges = interp1( Darwin_measures_centerline_distance, Darwin_measures_flowspeed, x_edges, 'linear', 'extrap');
plot(x_edges/1000, (measures_on_xedges - average_vel_estimate*(5/4))./measures_on_xedges, 'color', [0.7 0.7 0.7], 'linewidth', 1)
hold on
plot([x_edges(1) x_edges(end)]/1000, [0 0],'k--')
if (ii == 1)
 plot(x_edges/1000, (measures_on_xedges - average_vel_estimate*(5/4))./measures_on_xedges, 'k', 'linewidth', 2)
end  
ylabel({'Fractional';' velocity'; 'difference';' [obs-model]'}, 'fontweight', 'bold')
xlim([x_P(1) x_P(end)]/1000)
%ylim([-25 25])
ylim([-10 1])
xlabel('Distance along flowband (km)', 'fontweight', 'bold')


% Hatherton : (modern surface - calculated surface)
subplot('position', [0.2 0.33 0.8 0.15])
set(gca, 'fontsize', 16)
plot(x_P2/1000, (S_modern2-S_P2(1,:))./S_modern2, 'color', [0.7 0.7 0.7], 'linewidth', 1)
hold on
plot([x_P2(1) x_P2(end)]/1000, [0 0],'k--')
if (ii == 1)
 plot(x_P2/1000, (S_modern2-S_P2(1,:))./S_modern2, 'k', 'linewidth', 2)
end
ylabel({'Fractional';' elevation ';'difference ';'[obs-model]'}, 'fontweight', 'bold')
xlim([x_P2(1) x_P2(end)]/1000)
ylim([-0.2 0.2])
if (ii == 6)
text(60, -100,'Hatherton Glacier', 'fontangle', 'italic', 'fontsize', 16)
end
set(gca, 'xticklabel', [])

% Hatherton : (observed - calculated surface velocity)
subplot('position', [0.2 0.12 0.8 0.15])
set(gca, 'fontsize', 16)
measures_on_xedges2 = interp1( Hat_measures_centerline_distance, Hat_measures_flowspeed, x_edges2, 'linear', 'extrap');
plot(x_edges2/1000, (measures_on_xedges2 - average_vel_estimate2*(5/4))./measures_on_xedges2, 'color', [0.7 0.7 0.7], 'linewidth', 1)
hold on
plot([x_edges2(1) x_edges2(end)]/1000, [0 0],'k--')
if (ii == 1)
   plot(x_edges2/1000, (measures_on_xedges2 - average_vel_estimate2*(5/4))./measures_on_xedges2, 'k', 'linewidth', 2)
end
xlabel('Distance along flowband (km)', 'fontweight', 'bold')
ylabel({'Fractional';' velocity'; 'difference';' [obs-model]'}, 'fontweight', 'bold')
xlim([x_P2(1) x_P2(end)]/1000)
ylim([-10 1])
%set(gca, 'ytick', [-75 -50 -25 0], 'yticklabel',{'-75','-50','-25','0'})

end

export_fig figure1_match_to_modern.pdf -pdf -transparent
save figure1_match_to_modern.fig

