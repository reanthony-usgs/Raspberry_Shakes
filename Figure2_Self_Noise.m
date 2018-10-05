% Plots the PSDs and Self-Noise for the Raspberry Shakes 

clear all

% set grey color (between 0-1) larger numbers are closer to white

GC = 0.4;

% load in the data

load Figure_2_Self_Noise_Data/Epi_ACC.mat

load Figure_2_Self_Noise_Data/ShakeNoise.mat

load Figure_2_Self_Noise_Data/PSD_STS2.mat

load Figure_2_Self_Noise_Data/PSD_Compact.mat

load Figure_2_Self_Noise_Data/PSD_MEMS.mat

load lnm_tcmpvertical.mat

tcompact_self_noise = mediannm;

% Convert Adam's measurements of power to dB
dB_EpiSensor = 10*log10(Epimodelmnm(:,1));

EpiSensor = [Epimodelmnm(:,2),dB_EpiSensor]; 

% Remove Nans

nan_index = find(isnan(ShakesNoise(:,2)));
ShakesNoise(nan_index,:) = [];

nan_index = find(isnan(PSD_STS2(:,2)));
PSD_STS2(nan_index,:) = [];


nan_index = find(isnan(PSD_Compact(:,2)));
PSD_Compact(nan_index,:) = [];

nan_index = find(isnan(PSD_MEMS(:,2)));
PSD_MEMS(nan_index,:) = [];


%% Make the figure

Hz_lines = [0.02,0.05,0.1,0.2,0.5,1,3,5,10,20,30,50];

ticks = (Hz_lines);
HZ_label = (ticks);

%Make the Peterson curves
fs=250;
dlP=.05;
PSDTOL=15;
[LNMA,HNMA,lpd1,lpd2]=peterson_acc(dlP,fs);

%Smoothed Peterson curves for plotting
NMplotind=(0.001:dlP:10);
LNMAp=spline(lpd1,LNMA,NMplotind);
HNMAp=spline(lpd2,HNMA,NMplotind);

pd1 = 10.^(lpd1);
pd2 = 10.^(lpd2);

% Remove the "Fake" part of the LNM and HNM above 10 Hz

SI = find(pd1 == 0.1);


pd1 = pd1(SI:end);
pd2 = pd2(SI:end);
LNMA = LNMA(SI:end);
HNMA= HNMA(SI:end);

figure(8);clf 

H1 = semilogx(ShakesNoise(:,1),ShakesNoise(:,2),'color',[GC,GC,GC]);
hold on
H2 = semilogx(ShakesNoise(:,1),ShakesNoise(:,3),'m');
H5 = semilogx(PSD_Compact(:,1),PSD_Compact(:,2),'color',[0,0.5,0]);
H6 = semilogx(tcompact_self_noise(:,2), tcompact_self_noise(:,1), 'color', rgb('Olive'));
H7 = semilogx(PSD_STS2(:,1),PSD_STS2(:,2),'b');
H3 = semilogx(1./pd1,LNMA,'k:');
H4 = semilogx(1./pd2,HNMA,'k:');



%H5 = semilogx(log10(Periods),Calm_PSD,'b');
%H6 = semilogx(log10(Periods),Wind_PSD,'r');
%legend('High Discharge Median','Low Discharge Median')

set(gca,'FontSize',24)
xlim([0.05 40])
ylim([-190 -90])
set(gca,'ydir','normal')
lgd = legend('RS4D geophone Self-Noise', 'RS4D geophone PSD', 'Trillium Compact PSD', 'Trillium Compact Self-Noise', 'STS-2 PSD')
lgd.FontSize = 18


set(H4,'LineWidth',3.0);
set(H3,'LineWidth',3.0);
set(H1,'LineWidth',2.0);
set(H2,'LineWidth',2.0);
set(H5,'LineWidth',2.0);
set(H6,'LineWidth',2.0);
set(H7,'LineWidth',2.0);


set(gca,'xtick',ticks)
set(gca,'Xticklabel',HZ_label)
%set(gca,'xdir','reverse')

xlabel('Frequency (Hz)')
ylabel('dB (rel. 1 (m/s^2)^2/Hz)')

%% Make Figure for the Accelerometers 

figure(9);clf 

H2 = semilogx(PSD_MEMS(:,1),PSD_MEMS(:,2),'m');
hold on
H5 = semilogx(EpiSensor(:,1),EpiSensor(:,2),'color',[0,0.5,0]);
H6 = semilogx(PSD_STS2(:,1),PSD_STS2(:,2),'b');
H3 = semilogx(1./pd1,LNMA,'k:');
H4 = semilogx(1./pd2,HNMA,'k:');


%H5 = semilogx(log10(Periods),Calm_PSD,'b');
%H6 = semilogx(log10(Periods),Wind_PSD,'r');
%legend('High Discharge Median','Low Discharge Median')

set(gca,'FontSize',24)
xlim([0.05 40])
ylim([-190 -50])
set(gca,'ydir','normal')
lgd = legend('RS4D MEMS Self-Noise/PSD', 'EpiSensor Self-Noise', 'STS-2 PSD')
lgd.FontSize = 18


set(H4,'LineWidth',3.0);
set(H3,'LineWidth',3.0);
set(H2,'LineWidth',2.0);
set(H5,'LineWidth',2.0);
set(H6,'LineWidth',2.0);

set(gca,'xtick',ticks)
set(gca,'Xticklabel',HZ_label)
%set(gca,'xdir','reverse')

xlabel('Frequency (Hz)')
ylabel('dB (rel. 1 (m/s^2)^2/Hz)')








