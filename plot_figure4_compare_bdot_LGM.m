% Compare accumulation rate patterns and how influence SS surface at LGM


addpath('export_fig')



load LGM_steady_state_min1.mat
bdot_1 = b_dot_P;
bdot2_1 = b_dot_P2;
S_P_1  = S_P;
S_P2_1 = S_P2;

load min1_bdot2_LGM.mat
bdot_2 = b_dot_P;
bdot2_2 = b_dot_P2;
S_P_2  = S_P;
S_P2_2 = S_P2;

load min1_bdot3_LGM.mat
bdot_3 = b_dot_P;
bdot2_3 = b_dot_P2;
S_P_3  = S_P;
S_P2_3 = S_P2;

% load min1_bdot4_LGM.mat
% bdot_4 = b_dot_P;
% S_P_4  = S_P;
% S_P2_4 = S_P2;


figure(4)
set(gcf, 'Units', 'centimeters','position', [35 20 22 20])
subplot('position', [0.12 0.73 0.35 0.2])
plot(x_P/1000, bdot_1(1,:), 'k', 'linewidth', 2)
hold on
plot(x_P/1000, bdot_2(1,:), 'c', 'linewidth', 2)
plot(x_P/1000, bdot_3(1,:), 'b', 'linewidth', 2)
ylabel({'Accumulation';' rate (m/yr)'}, 'fontweight', 'demi')
set(gca, 'fontsize', 14)
xlim([x_P(1) x_P(end)]/1000)
%xlabel('Distance along flowband (km)')
title ('Darwin Glacier')
legend('Modern', 'LGM estimate 1', 'LGM estimate 2', 'location', 'southeast', 'fontsize', 14)

subplot('position', [0.62 0.73 0.35 0.2])
plot(x_P2/1000, bdot2_1, 'k', 'linewidth', 2)
hold on
plot(x_P2/1000, bdot2_2, 'c', 'linewidth', 2)
plot(x_P2/1000, bdot2_3, 'b', 'linewidth', 2)
ylabel({'Accumulation';' rate (m/yr)'}, 'fontweight', 'demi')
set(gca, 'fontsize', 14)
xlim([x_P2(1) x_P2(end)]/1000)
%xlabel('Distance along flowband (km)')
title('Hatherton Glacier')



subplot('position', [0.12 0.42 0.35 0.24])
plot(x_P/1000, S_P_1, 'k')
hold on
plot(x_P/1000, S_P_2, 'c')
plot(x_P/1000, S_P_3, 'b')
ylabel({'Elevation (m)'}, 'fontweight', 'demi')
set(gca, 'fontsize', 14)
xlim([x_P(1) x_P(end)]/1000)
ylim([800 2200])

subplot('position', [0.62 0.42 0.35 0.24])
plot(x_P2/1000, S_P2_1, 'k')
hold on
plot(x_P2/1000, S_P2_2, 'c')
plot(x_P2/1000, S_P2_3, 'b')
ylabel({'Elevation (m)'}, 'fontweight', 'demi')
set(gca, 'fontsize', 14)
xlim([x_P2(1) x_P2(end)]/1000)



subplot('position', [0.12 0.12 0.35 0.24])
plot(x_P/1000, S_P_2(1,:) - S_P_1(1,:), 'c')
hold on
plot(x_P/1000, S_P_3(1,:)-S_P_1(1,:),'b')
ylabel({'Elevation';'difference (m)'}, 'fontweight', 'demi')
set(gca, 'fontsize', 14)
xlim([x_P(1) x_P(end)]/1000)
xlabel('Distance along flowband (km)')

subplot('position', [0.62 0.12 0.35 0.24])
plot(x_P2/1000, S_P2_2(1,:) - S_P2_1(1,:), 'c')
hold on
plot(x_P2/1000, S_P2_3(1,:)-S_P2_1(1,:),'b')
ylabel({'Elevation';'difference (m)'}, 'fontweight', 'demi')
set(gca, 'fontsize', 14)
xlim([x_P2(1) x_P2(end)]/1000)
xlabel('Distance along flowband (km)')


export_fig figure4_compare_bdot.pdf -pdf -transparent


% LGM_steady_state_min2.mat
% min2_bdot2_LGM.mat
% min2_bdot3_LGM.mat
% min2_bdot4_LGM.mat

