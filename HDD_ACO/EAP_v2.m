% HDD Test Cases Code - > Franklin Wash

close all
clear

folder = '/Users/lfp3/Dropbox (GaTech)/GT/Fall_19/HDD-ACO/Case Studies to be evaluated/AEP/';

ACO_1 = readmatrix(fullfile(folder,'ACO.csv'));
ACO_2 = readmatrix(fullfile(folder,'ACO_2.csv'));
AD_Alignment = readmatrix(fullfile(folder,'Alignment.csv')); % Assumes Actual design was there 
Alignment = AD_Alignment;
Profile = readmatrix(fullfile(folder,'Profile.csv'));
Profile(end,:)=[];

%% Plot ACO evolution

figure(10)
%subplot(2,1,1)
hold on

nG = size(ACO_1,1);
cf=autumn(nG+10);
colors=flipud(cf(1:nG,:));

% ACO_1
ACO = ACO_1;

for g=1:nG
    
    HDD_Entry = [ACO(g,8) interp1(Profile(:,1),Profile(:,2),ACO(g,8))];
    Alignment(1,2)=ACO(g,5); %alpha
    Alignment(7,2)=ACO(g,6); %beta
    Alignment(1,1)=ACO(g,9); %LI
    Alignment(7,1)=ACO(g,10); %LIX
    %Alignment(3,1)=ACO(g,11); %LM
    
    [a]=plotHDD(Alignment,HDD_Entry,colors(g,:));
end
a=plotHDD(Alignment,HDD_Entry,[60 0 15]./255);

% ACO_1
ACO = ACO_2;

for g=1:nG
    
    HDD_Entry = [ACO(g,8) interp1(Profile(:,1),Profile(:,2),ACO(g,8))];
    Alignment(1,2)=ACO(g,5); %alpha
    Alignment(7,2)=ACO(g,6); %beta
    Alignment(1,1)=ACO(g,9); %LI
    Alignment(7,1)=ACO(g,10); %LIX
    %Alignment(3,1)=ACO(g,11); %LM
    
    [a]=plotHDD(Alignment,HDD_Entry,colors(g,:));
end
a=plotHDD(Alignment,HDD_Entry,[60 0 15]./255);

%% Plot Actual design

figure(10)
hold on
plot(Profile(:,1),Profile(:,2),'k-')

Entry_station = AD_Alignment(1,16)+50;
HDD_Entry = [Entry_station interp1(Profile(:,1),Profile(:,2),Entry_station)];
plotHDD(AD_Alignment,HDD_Entry,[0 38 154]./255);

ylim([-110 20])
xlim([0 2050])
xlabel('Stations')
ylabel('Elevation [ft]')

ax = gca;
ax.XRuler.Exponent = 0;

%% Plot Length evolution 

GS = 0:nG-1;

% ACO_1
ACO = ACO_1;

figure(15);

hold on

%L(4)=plot([0 30],[1960 1960],'Color',[0 38 154]./255,'LineStyle','--'); % Actual Design
plot([0 30],[1960 1960],'Color',[0 38 154]./255,'LineStyle','--'); % Actual Design


h11=plot(GS,1000./ACO(:,1),'Color',[0 50 82]./255); % Best
h12=plot(GS,1000./ACO(:,2),'Color',[177 83 15]./255,'LineStyle','-.'); % Mean
h13=plot(GS,1000./ACO(:,3),'Color',[250 129 22]./255,'LineStyle',':'); % Q25
plot(GS,1000./ACO(:,4),'Color',[250 129 22]./255,'LineStyle',':'); % Q75

% ACO_2
ACO = ACO_2;

figure(15)
hold on
h2(1)=plot(GS,1000./ACO(:,1),'Color',[0.6350 0.0780 0.1840]); % Best
h2(2)=plot(GS,1000./ACO(:,2),'Color',[0 50 0]./255,'LineStyle','-.'); % Mean
h2(3)=plot(GS,1000./ACO(:,3),'Color',[0.4660 0.6740 0.1880],'LineStyle',':'); % Q25
plot(GS,1000./ACO(:,4),'Color',[0.4660 0.6740 0.1880],'LineStyle',':'); % Q75

legend_entries={'Mean Length','Minimum Length','1st and 3rd Quartiles'};
%legend(h2(1,2,4),legend_entries);

% a=axes('position',get(gca,'position'),'visible','off');
% legend(a,h2,legend_entries,'Location','EastOutside');

L(1) = plot(nan, nan,'Color','k');
L(2) = plot(nan, nan ,'LineStyle','-.','Color','k');
L(3) = plot(nan, nan,'LineStyle',':','Color','k');
legend(L,legend_entries)

%ylim([2252 2270])
title('Total Length Evolution')
xlabel('Generation')
ylabel('Total Drill length [ft]')

%% Convergence Evolution

figure(18)
subplot(1,2,1)
hold on
ml = min(ACO_2(:,1));
plot(GS,abs(ACO_2(:,4)-ACO_2(:,3))./ml,'Color',[0 50 0]./255) % IQR
plot(GS,abs(ACO_2(:,2)-ACO_2(:,1))./ml,'Color',[0.6350 0.0780 0.1840]) % Mean-Min
xlabel('Generation')
ylabel('Normalized Deviation [%]')
title('ACO - 50 ft')
ytickformat('percentage')
ylim([1e-4 10])
set(gca, 'YScale', 'log')

subplot(1,2,2)
hold on
ml = min(ACO_1(:,1));
plot(GS,abs(ACO_1(:,4)-ACO_1(:,3)).*(100/ml),'Color',[0 50 0]./255) % IQR
plot(GS,abs(ACO_1(:,2)-ACO_1(:,1)).*(100/ml),'Color',[0.6350 0.0780 0.1840]) % Mean-Min
xlabel('Generation')
ylabel('Normalized Deviation')
title('ACO - 150 ft')
ytickformat('percentage')
ylim([1e-4 10])
set(gca, 'YScale', 'log')

legend('Interquartile Range','Mean-Min Difference')


figure(16)
Names={'Entry Angle','Exit Angle','Bottom path Elevation'};
ylabels={'Angle [deg]','Angle [deg]','Elevation [ft]'};

for pl=1:3
    subplot(3,2,((pl-1)*2+1))
    hold on
    plot(GS,ACO_2(:,4+pl),'Color',[0.6350 0.0780 0.1840]) % Best
    plot(GS,ACO_2(:,11+pl),'Color',[0 50 0]./255) % Mean
    %title(Names{pl})   
    ylabel(ylabels{pl})
    xlabel('Generations')
    
    subplot(3,2,(pl-1)*2+2)
    hold on
    plot(GS,ACO_1(:,4+pl),'Color',[0.6350 0.0780 0.1840]) % Best
    plot(GS,ACO_1(:,11+pl),'Color',[0 50 0]./255) % Mean
    %title(Names{pl})   
    ylabel(ylabels{pl})
    xlabel('Generations')
    
end

legend('Minimum Length','Mean Length')