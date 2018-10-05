% Initial statistical Processing of the Oklahoma Trillium Compacts 
clear all

load OK_TCompact_SNR.mat
                  
% Get the Difference
Diff_Mag = SNR_T_Compacts_mat(:,1) - SNR_T_Compacts_mat(:,2);

OK_TCompact = [Diff_Mag, SNR_T_Compacts_mat(:,2:end)];


% set range parameters 
min_range =20;
max_range = 40;
Big_Event_Threshold =3;
Range_Thres = 167;


% remove all data with events ouside of 1 degree from the station

too_far = find(OK_TCompact(:,3) > Range_Thres);
OK_TCompact(too_far,:) = [];



% data format (misfit, NEIC ML, EpiCenter Distance)

% Find the mean misfit and correct since we are using vertical components
% (Empirical Correction) 

Mean_misfit = mean(OK_TCompact(:,1))

OK_TCompact(:,1) = OK_TCompact(:,1);

std(OK_TCompact(:,1))


Big_Events_index = find(OK_TCompact(:,2) >= Big_Event_Threshold);

Big_Events = OK_TCompact(Big_Events_index,:);

display('Big Event Mean and STD')

Mean_misfit = mean(Big_Events(:,1))

std(Big_Events(:,1))






% make a histogram of misfits 

figure(1); clf 
histogram(OK_TCompact(:,1),40, 'BinLimits', [-2,2], 'Normalization', 'probability')
hold on
histogram(Big_Events(:,1), 40, 'BinLimits', [-2,2], 'Normalization', 'probability')
xlabel('Misfit (Station M_L-NEIC M_L)')
ylabel('Probability')

h= legend('All Events', 'M$_{L}$ $\geq$ 3.0')
set(h,'Interpreter','latex')
ax = gca;
set(gca,'FontSize',28)
%c = ax.Color;
ax.LineWidth = 3;
h= legend('All Events', 'M$_{L}$ $\geq$ 3.0')
set(h,'Interpreter','latex')

%% Scatter plot of misfit with distance and magnitude 



figure(2);clf
set(gca,'FontSize',30)

scatter(OK_TCompact(:,3),OK_TCompact(:,2),100,OK_TCompact(:,1),'filled')
colormap(redblue);
c=colorbar
xlabel('Station Distance from EpiCenter')
ylabel('NEIC ML')
caxis([-1 1])
%ylabel(c,'Station ML Misfit')

ax = gca;
%c = ax.Color;
ax.LineWidth = 3;


%% Ok let's plot SNR vs Magnitude for a range

In_Range_Index = find(OK_TCompact(:,3) >= min_range & OK_TCompact(:,3) <= max_range);

Range_Events = OK_TCompact(In_Range_Index,:);

SNR = Range_Events(:,8);

Edges = (1.95:0.1:4.05);

Edge_Plot = (2:0.1:4.0);

[N,Edges2,bins] = histcounts(Range_Events(:,2),Edges);

for kk = 1:length(Edges)-1
    if N(kk) >= 10
        power_median(kk) = nanmedian(SNR(bins==kk));
    else
        power_median(kk) = NaN;
    end

end





figure(3);clf
set(gca,'FontSize',30)
plot(Range_Events(:,2),Range_Events(:,8),'x')
hold on
plot(Edge_Plot,power_median,'o-','linewidth',3)
xlabel('NEIC Magnitude')
ylabel('SNR')

ylim([0 12])

ax = gca;
set(gca,'FontSize',20)
%c = ax.Color;
ax.LineWidth = 3;


USGS_SNR_20_40 = [Edge_Plot', power_median'];



%% Also SNR vs MISfit 

%figure(4);clf
%set(gca,'FontSize',30)
%plot(Range_Events(:,8),Range_Events(:,1),'x')
%xlabel('SNR')
%ylabel('Misfit (Station M_L-NEIC M_L)')
%ylim([-2 5])

%ax = gca;
%set(gca,'FontSize',20)
%c = ax.Color;
%ax.LineWidth = 3;



