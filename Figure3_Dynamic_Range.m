% Plots the PSDs and Self-Noise for the Raspberry Shakes 

clear all

% set grey color (between 0-1) larger numbers are closer to white

GC = 0.4;

% load in the data

% All files are in format of [frequency, amplitude} 

%load Figure_2_Self_Noise_Data/Epi_ACC.mat

load Figure_2_Self_Noise_Data/ShakeNoise.mat

load Figure_2_Self_Noise_Data/PSD_STS2.mat

load Figure_2_Self_Noise_Data/PSD_Compact.mat

load Figure_2_Self_Noise_Data/PSD_MEMS.mat

% Convert Adam's measurements to dB
%dB_EpiSensor = 20*log10(Epimodelmnm(:,1));

% make Adam's file look like the others
%EpiSensor = [Epimodelmnm(:,2),dB_EpiSensor]; 

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

%Convert from dB back to ground acceleration
ShakesNoise(:,2) = 10.^(ShakesNoise(:,2)./10);
% Now in (m/s/s)^2/Hz
ShakesNoise(:,2) = sqrt(ShakesNoise(:,2).*ShakesNoise(:,1).*(2^(1/2)-2^(-1/2)));
% Now in one-octave integrated bands


H1 = loglog(PSD_MEMS(:,1),ShakesNoise(:,2),'color',[GC,GC,GC]);




hold on

PSD_MEMS(:,2) = sqrt(10.^(PSD_MEMS(:,2)./10).*(PSD_MEMS(:,1))*(2^(1/2)-2^(-1/2)));

%H5 = loglog(PSD_MEMS(:,1),PSD_MEMS(:,2),'k');


% Should be in amplitude of sine so add 


% You want an amplitude rms so sqrt(2) needs to be multiplied.
%H2 = loglog(PSD_MEMS(:,1),0.707*PSD_MEMS(:,1)*2*pi*0.022,'--','color',[GC,GC,GC]);
%H6 = loglog(PSD_MEMS(:,1),0.707*ones(length(PSD_MEMS),1)*20,'k--');


%LNMA = integrateme(LNMA, 1./pd1,
%HNMA = integrateme(HNMA, 1./pd2, 


%LNMA = sqrt((10.^(LNMA./10)).*(1./pd1));
%HNMA = sqrt((10.^(HNMA./10)).*(1./pd2));

LNMA = sqrt((10.^(LNMA./10)).*(1./pd1)*(2^(1/2)-2^(-1/2)));
HNMA = sqrt((10.^(HNMA./10)).*(1./pd2)*(2^(1/2)-2^(-1/2)));



set(gca,'FontSize',24)
xlim([0.05 40])
ylim([10^(-10) 10^3])
set(gca,'ydir','normal')
%lgd = legend('RS-4D geophone Self-Noise', 'RS-4D MEMS Self-Noise', 'RS-4D Geophone Clip Level', 'RS-4D MEMS Clip Level')
%lgd.FontSize = 18





%set(gca,'xdir','reverse')

xlabel('Frequency (Hz)')
ylabel('Octave Wide Bandpassed Acceleration (m/s^2)')


% Let's try this as a colored area map

One_M = ones(355,1)';


X = PSD_MEMS(1:355,1);

Y = [ShakesNoise(1:355,2)'; PSD_MEMS(1:355,2)'; sqrt(0.707*PSD_MEMS(1:355,1)'*(2*pi*0.022)^2); sqrt(.707*9.8*9.8.*One_M); 0.707*One_M*10000000]

Y=Y';
a = area(X,Y)
set (gca, 'Xscale', 'log')
alpha(0.3)

a(1).FaceColor = 'r'
a(2).FaceColor = 'c'
a(3).FaceColor = 'g'
a(4).FaceColor = 'y'
a(5).FaceColor = 'r'


H1 = loglog(PSD_MEMS(:,1),ShakesNoise(:,2),'k');
H5 = loglog(PSD_MEMS(:,1),PSD_MEMS(:,2),'color', rgb('Purple'));
H2 = loglog(PSD_MEMS(:,1),sqrt(0.707*PSD_MEMS(:,1)*(2*pi*0.022).^2),'color',[GC,GC,GC]);
H6 = loglog(PSD_MEMS(:,1),sqrt(0.707*ones(length(PSD_MEMS),1)*9.8*9.8),'b-');

H3 = loglog(1./pd1,LNMA,'k:');
H4 = loglog(1./pd2,HNMA,'k:');


set(H1,'LineWidth',4.0);
set(H2,'LineWidth',4.0);
set(H5,'LineWidth',4.0);
set(H6,'LineWidth',4.0);


set(H4,'LineWidth',3.0);
set(H3,'LineWidth',3.0);

ax = gca;
%c = ax.Color;
ax.LineWidth = 3;
set(gca,'xtick',ticks)
set(gca,'Xticklabel',HZ_label)
set(gca,'Layer','top')

lgd = legend([H1,H2,H5,H6],{'RS-4D Geophone Self-Noise', 'RS-4D Geophone Clip Level'  'RS-4D MEMS Self-Noise' , 'RS-4D MEMS Clip Level'})
lgd.FontSize = 18



