% HDD Test Cases Code - > Franklin Wash

%close all

folder = '/Users/lfp3/Dropbox (GaTech)/GT/Fall_19/HDD-ACO/Case Studies to be evaluated/AEP/';

ACO = readmatrix(fullfile(folder,'ACO_2.csv'));
AD_Alignment = readmatrix(fullfile(folder,'Alignment.csv')); % Assumes Actual design was there 
Alignment = AD_Alignment;
Profile = readmatrix(fullfile(folder,'Profile.csv'));
Profile(end,:)=[];

%% Plot ACO evolution

figure(10)
%subplot(2,1,1)
hold on

nG = size(ACO,1);
cf=autumn(nG+10);
colors=flipud(cf(1:nG,:));

for g=1:nG
    
    HDD_Entry = [ACO(g,8) interp1(Profile(:,1),Profile(:,2),ACO(g,8))];
    Alignment(1,2)=ACO(g,5); %alpha
    Alignment(7,2)=ACO(g,6); %beta
    Alignment(1,1)=ACO(g,9); %LI
    Alignment(7,1)=ACO(g,10); %LIX
    %Alignment(3,1)=ACO(g,11); %LM
    
    [a]=plotHDD(Alignment,HDD_Entry,colors(g,:));
end

%% Plot Actual design

figure(10)
%subplot(2,1,1)
hold on
plot(Profile(:,1),Profile(:,2),'k-')
plotHDD(Alignment,HDD_Entry,[60 0 15]./255);

Entry_station = AD_Alignment(1,16)+50;
HDD_Entry = [Entry_station interp1(Profile(:,1),Profile(:,2),Entry_station)];

plotHDD(AD_Alignment,HDD_Entry,[0 38 154]./255);

ylim([-110 20])
xlim([0 2050])
xlabel('Stations')
ylabel('Elevation [ft]')

ax = gca;
ax.XRuler.Exponent = 0;

%% Plot Convergence evolution 

GS = 0:nG-1;

figure(15)
subplot(2,1,1);
hold on
plot(GS,1000./ACO(:,1),'Color',[0.6350 0.0780 0.1840]) % Best
plot(GS,1000./ACO(:,2),'Color',[0 50 0]./255) % Mean
plot(GS,1000./ACO(:,3),'Color',[0.4660 0.6740 0.1880],'LineStyle',':') % Q25
plot(GS,1000./ACO(:,4),'Color',[0.4660 0.6740 0.1880],'LineStyle',':') % Q75
title('Total Length Evolution')

plot([0 30],[1600 1600],'Color',[0.8500 0.3250 0.0980]) % Geometric Minimum
plot([0 30],[1960 1960],'Color',[0 38 154]./255) % Actual Design

%ylim([2252 2270])
xlabel('Generation')
ylabel('Total Drill length [ft]')
legend_entries={'Engineered Design','Geometric Minimum','1st and 4th  Quartiles','Mean Length','Minimum Length'};

h=get(gca,'Children'); % grab all the axes handles at once
legend(h([1 2 3 5 6]),legend_entries) 

subplot(2,1,2)
hold on
plot(GS,abs(ACO(:,4)-ACO(:,3)).*(1000/2254),'Color',[0 50 0]./255) % IQR
plot(GS,abs(ACO(:,2)-ACO(:,1)).*(1000/2254),'Color',[0.6350 0.0780 0.1840]) % Mean-Min
legend('Interquartile Range','Mean-Min Difference')
xlabel('Generation')
ylabel('Normalized Deviation')

figure(16)
Names={'Entry Angle','Exit Angle','Bottom path Elevation'};
ylabels={'Angle [deg]','Angle [deg]','Elevation [ft]'};

for pl=1:3
    subplot(3,1,pl)
    hold on
    plot(GS,ACO(:,4+pl),'Color',[0.6350 0.0780 0.1840]) % Best
    plot(GS,ACO(:,11+pl),'Color',[0 50 0]./255) % Mean
    title(Names{pl})   
    ylabel(ylabels{pl})
    xlabel('Generations')
    
end

legend('Minimum Length','Mean Length')