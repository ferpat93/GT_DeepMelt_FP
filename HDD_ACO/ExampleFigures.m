% HDD Test Cases Code - > Franklin Wash

close all

folder = '/Users/lfp3/Dropbox (GaTech)/GT/Fall_19/HDD-ACO/';


%% Case 2

Case1_1 = readmatrix(fullfile('Case1_1.csv'));
Case1_2 = readmatrix(fullfile('Case1_2.csv'));
Profile = readmatrix(fullfile(folder,'ExampleProfile.csv'));

figure(11)
subplot(2,1,2)
hold on
plot(Profile(:,1),Profile(:,2),'k:') % Profile

Entry_station = Case1_1(1,16);
HDD_Entry = [Entry_station interp1(Profile(:,1),Profile(:,2),Entry_station)];

plotHDD(Case1_1,HDD_Entry,[1 0 0])

Entry_station = Case1_2(1,16);
HDD_Entry = [Entry_station interp1(Profile(:,1),Profile(:,2),Entry_station)];

plotHDD(Case1_2,HDD_Entry,[0 0 1])

legend("Ground Surface","Drill Path 1","Drill Path 2")
title('Case 2: Fixed Intermediate Section')

ylim([0 110])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

%% Case 1

Case2_1 = readmatrix(fullfile('Case2_1.csv'));
Case2_2 = readmatrix(fullfile('Case2_2.csv'));
Profile = readmatrix(fullfile(folder,'ExampleProfile.csv'));

Entry_station = 100;
HDD_Entry = [Entry_station interp1(Profile(:,1),Profile(:,2),Entry_station)];

figure(11)
subplot(2,1,1)
hold on
plot(Profile(:,1),Profile(:,2),'k:') % Profile
plotHDD(Case2_1,HDD_Entry,[1 0 0])
plotHDD(Case2_2,HDD_Entry,[0 0 1])

legend("Ground Surface","Drill Path 1","Drill Path 2")
title('Case 1: Fixed Bounds (Entry/Exit points)')

ylim([0 110])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])