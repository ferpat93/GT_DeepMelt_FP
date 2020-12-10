% HDD Test Cases Code - > FixedBounds
close all

folder =  '/Users/lfp3/Dropbox (GaTech)/GT/Previous Semesters/Fall_19/HDD-ACO/Case Studies to be evaluated/FixedBounds/';

gens = [0 5 15 30 40 50];
ACO = readmatrix(fullfile(folder,'ACO.xlsx'));

Profile = [(0:10:2000)' 100*ones(201,1)];

Short = readmatrix(fullfile(folder,'MinMax.xlsx'),'Sheet','short');
Long = readmatrix(fullfile(folder,'MinMax.xlsx'),'Sheet','long');

for gi=1:numel(gens)
    Ledger(:,:,gi) = readmatrix(fullfile(folder,'LedgerOverGen.xlsx'),'Sheet',num2str(gens(gi)));
end

Entry_station = Short(1,16);
HDD_Entry = [Entry_station interp1(Profile(:,1),Profile(:,2),Entry_station)];

%% Plot ACO evolution

figure(10)
%subplot(1,2,1)
hold on

nG = size(ACO,1);
cf=autumn(nG+10);
colors=flipud(cf(1:nG,:));
Alignment = Short;

for g=1:nG
    
    HDD_Entry = [ACO(g,8) interp1(Profile(:,1),Profile(:,2),ACO(g,8))];
    Alignment(1,2)=ACO(g,5); %alpha
    Alignment(5,2)=ACO(g,6); %beta
    Alignment(1,1)=ACO(g,9); %LI
    Alignment(5,1)=ACO(g,10); %LIX
    Alignment(3,1)=ACO(g,11); %LM
    
    a=plotHDD(Alignment,HDD_Entry,colors(g,:));
end

%% Plot Actual design
figure(10)
%subplot(1,2,1)
hold on
plot(Profile(:,1),Profile(:,2),'k-')
plotHDD(Short,HDD_Entry,[60 0 15]./255);
plotHDD(Long,HDD_Entry,[15 0 60]./255);

xlabel('Stations')
ylabel('Elevation [ft]')
xlim([0 2000])
ylim([-200 120])

% Entry_station = AD_Alignment(1,16)+0;
% HDD_Entry = [Entry_station interp1(Profile(:,1),Profile(:,2),Entry_station)];
% 
% plotHDD(AD_Alignment,HDD_Entry,[0 38 154]./255);
% 
% ylim([-10 100])
% 
% 
% ax = gca;
% ax.XRuler.Exponent = 0;

%% PLOT Alpha distribution
bins = (10:1:20)';
counts = zeros(numel(bins)-1,numel(gens));
for gi=1:numel(gens)
    [N,~] = histcounts(Ledger(:,3,gi),bins);
    counts(:,gi)=N';
end

xvalues = string(gens);
yvalues = string(bins(1:end-1));
% , 'Colormap', [1 1 1 ; flipud(summer)]
figure;
h = heatmap(xvalues,yvalues,counts);

h.Title = 'Distribution of entry angles over generations';
h.XLabel = 'Generations';
h.YLabel = 'Alpha [deg]';
%h.ColorScaling = 'scaledcolumns';
%% Plot Convergence evolution 

GS = 0:nG-1;
minL = Short(end,14);
maxL = Long(end,14);
rangeL = maxL-minL;

figure(15)
hold on
plot(GS,((1000./ACO(:,1))-minL)/rangeL,'Color',[0.6350 0.0780 0.1840]) % Best
plot(GS,((1000./ACO(:,2))-minL)/rangeL,'Color',[0 50 0]./255) % Mean
%plot(GS,((1000./ACO(:,3))-minL)/rangeL,'Color',[0.4660 0.6740 0.1880],'LineStyle',':') % Q25
%plot(GS,((1000./ACO(:,4))-minL)/rangeL,'Color',[0.4660 0.6740 0.1880],'LineStyle',':') % Q75
title('Normalized Length Evolution - Fixed Entry/Exit Points')

plot([0 50],[0 0],'Color',[0.8500 0.3250 0.0980]) % Geometric Minimum
plot([0 50],[1 1],'Color',[0 38 154]./255) % Actual Design

%ylim([0 1])
xlabel('Generation')
ylabel('Normalized Drill length [-]')
legend_entries={'Average Ledger Solution','Best ACO Solution'};
% 'Ledger 1st and 3rd  Quartiles ','Average Ledger Solution','Maximum Length','Minimum Length',
h=get(gca,'Children'); % grab all the axes handles at once
legend(h([3 4]),legend_entries) 
%set(gca, 'YScale', 'log')
grid on

