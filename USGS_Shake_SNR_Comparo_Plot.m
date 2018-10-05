% Comparison Plot of Shakes vs. USGS Afterschock Deployment SNRs

clear all

% load in the data

load USGS_SNR_20_40km.mat
load USGS_SNR_80_100km.mat
load Shakes_SNR_20_40km.mat
load Shakes_SNR_80_100km.mat

figure(7);clf
set(gca,'FontSize',30)
plot(USGS_SNR_20_40(:,1),USGS_SNR_20_40(:,2),'^-','linewidth',3, 'color',[0,0.5,0], 'MarkerEdgeColor','k','MarkerFaceColor',[0,0.5,0], 'MarkerSize',10)
hold on
plot(USGS_SNR_80_100(:,1),USGS_SNR_80_100(:,2),'^-','linewidth',2, 'color',[0,0.5,0],'MarkerEdgeColor','k','MarkerFaceColor',[0,0.5,0], 'MarkerSize',7)
plot(Shakes_SNR_20_40(:,1),Shakes_SNR_20_40(:,2),'m^-','linewidth',3,'MarkerEdgeColor','k','MarkerFaceColor','m', 'MarkerSize',10)
plot(Shakes_SNR_80_100(:,1),Shakes_SNR_80_100(:,2),'m^-','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','m', 'MarkerSize',7)
plot(Shakes_SNR_80_100(:,1), 2*ones(1,length(Shakes_SNR_80_100)),'k','linewidth',1)
xlabel('NEIC Magnitude (M_L)')
ylabel('Signal-to-Noise Ratio')
legend('USGS 20 to 40 km', 'USGS 80 to 100 km', 'Raspberry Shakes 20 to 40 km', 'Raspberry Shakes 80 to 100 km') 

ylim([0 12])

ax = gca;
set(gca,'FontSize',24)
%c = ax.Color;
ax.LineWidth = 3;



