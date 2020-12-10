close all
% 
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.065 0.065], [0.2 0.05], [0.075 0.04]); % subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
% if ~make_it_tight,  clear subplot;  end

%% ~ PRE - FUSION ~ %%

% Constants
Ao_mm = (13/20)^2*pi(); %mm2
cmap = [27,158,119; 217,95,2; 117,112,179]./255; % Colorblind safe colormap
prob=0.95; % Probability for confidence interval
A_int = 1:0.5:4; % Area Intervals for part 1 (relative to initial area)
DishTreatments = {'Control','Control','Glucose','Glucose','NaCl','NaCl',...
    'Control','Control','Glucose','Glucose','NaCl','NaCl',...
    'Glucose','Glucose','Glucose','Glucose',...
    'NaCl','NaCl','NaCl','NaCl'};

Treatments = {'Control','Glucose','NaCl'};
LegendTreatments = {'Neutral','Nutritive','Adverse'};

DataFolder = 'D:\Slime_Mold_Network\Data\Fusion';
ResultsFolder = 'D:\Slime_Mold_Network\Results\Fusion';

DataFolder = '/Users/lfp3/Dropbox (GaTech)/GT/Spring-20/SMN/Data';
ResultsFolder = '/Users/lfp3/Dropbox (GaTech)/GT/Spring-20/SMN/Data';

LedgerName = 'Ledger.xlsx';
Ledger = readtable(fullfile(DataFolder,LedgerName));

% Gather Data
gather = false;

if gather
    Folders=GetSubfolders(DataFolder);

    varTypes = [repelem({'double'},16)  {'string'}];
    TableVariables = {'cell','An','t','Ao','At','pC','Ae','Aws','nwsr','Wbins','Lbins','Wpar','Lpar','nNodes_nEdges','Le_p','nDish','folder'};
    sz = [0 numel(TableVariables)];
    T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',TableVariables);

    for f=1:numel(Folders)

        if isfile(fullfile(DataFolder,Folders{f},'PreFusion_indexes.mat'))
            Data_fi = extractfield(load(fullfile(DataFolder,Folders{f},'PreFusion_indexes.mat')),'Data_f');
            T = vertcat(T,unique(Data_fi{1},'rows'));
        end
    end

    T = addvars(T,DishTreatments(T.nDish)','NewVariableNames','Treatment');
    T(T.t==inf,:)=[]; % Erase null entries
    T.t=T.t./60;
    T(and(T.t==0,T.An>1),:)=[];
    T.Ae = Ao_mm.*T.Ae./T.Ao;
    T.At = Ao_mm.*T.At./T.Ao;
    T.Aws = Ao_mm.*T.Aws./T.Ao;
    T.Ao = Ao_mm.*ones(size(T.Ao));
    T.meanVein = 10.*T.An./T.Le_p(:,2); T(T.meanVein==inf,:)
    
    writetable(T,fullfile(DataFolder,'Data_BF.xlsx'),'Sheet',1,'Range','A1')
else
    T = readtable(fullfile(DataFolder,'Data_BF.xlsx'));
    T(and(T.t==0,T.An>1),:)=[];
    T.meanVein(T.meanVein==65535)=nan;
end

%% Reviewer Stuff %%


T.Tort = ((T.Le_p_2./T.Le_p_1)-1).*100;

for t=1:3 % Gather Data
    Tt = T(strcmp(T.Treatment,Treatments{t}),:);    
    x = Tt.t; y = Tt.Tort;
    
    figure(111); hold on
    scatter(Tt.t,Tt.Tort)
    
    figure(112); hold on
    scatter(Tt.An,Tt.Tort)
end

figure(111); legend(Treatments)
ylim([0 50]); xlim([0.1 inf]);
xlabel('Time [h]'); ylabel('Tortuosity [%]')
figure(112); legend(Treatments)
ylim([0 50]); xlim([1 inf])
xlabel('Normalized Area [-]'); ylabel('Tortuosity [%]')

%% -- START PLOTTING -- %%

% 1) Area vs time

Data = zeros(3,numel(A_int),numel(Treatments));

for t=1:3 % Gather Data
    Tt = T(strcmp(T.Treatment,Treatments{t}),:);
    
    for ati=1:numel(A_int)
        Tta = Tt(Tt.An==A_int(ati),:);
        Data(:,ati,t)=[mean(Tta.t) CI(Tta.t,prob)] ; 
        ppp=1;
    end
    
end

% figure;
% plotTreatmentsTime(A_int,Data,'Time Elapsed to target areas',...
%     'Normalized Area [-]','Time [h]',LegendTreatments)

% saveas(gcf,'BF_AreaVsTime.png')
% print('BF_AreaVsTime','-depsc');
% saveas(gcf,'BF_AreaVsTime.fig')

% figure; hold on
% y = A_int;
% yr = [y fliplr(y)];
% for t=1:3 % Plot Data
%     x = Data(1,:,t);
%     xr = [x+Data(2,:,t), fliplr(x+Data(3,:,t))];
%     f=fill(xr,yr,cmap(t,:));
%     set(f,'facealpha',.2,'LineStyle','none')
%     p(t) = plot(x,y, 'Color',cmap(t,:), 'LineWidth', 2);
% end
% xlabel('Time [h]')
% ylabel('Normalized Area [-]')
% title()
% legend([p(1) p(2) p(3)],Treatments)



% 2) SM Density

% Data = zeros(3,numel(A_int),numel(Treatments));
% 
% for t=1:3 % Gather Data
%     Tt = T(strcmp(T.Treatment,Treatments{t}),:);
%     
%     for ati=1:numel(A_int)
%         Tta = Tt(Tt.An==A_int(ati),:);
%         Dens = Tta.At./Tta.Ae;
%         Data(:,ati,t)=[mean(Dens) CI(Dens,prob)] ; 
%     end
%     
% end
% 
% plotTreatmentsTime(A_int,Data,'Slime mold Solidity','Normalized Area [-]','Slime Mold Solidity [%]',Treatments)
% saveas(gcf,'BF_Solidity.png')
% saveas(gcf,'BF_Solidity.fig')

% 3) SM Whitespace Area mean

Data = zeros(3,numel(A_int),numel(Treatments));
DataRWS = zeros(3,numel(A_int),numel(Treatments));

MeanEnclosed = 100.*T.Aws./T.nwsr; MeanEnclosed(isinf(MeanEnclosed)) = nan;
T.MeanEnclosed = MeanEnclosed;
T.RatioEmpty = T.Aws./T.At;

for t=1:3 % Gather Data
    Tt = T(strcmp(T.Treatment,Treatments{t}),:);
    
    for ati=1:numel(A_int)
        Tta = Tt(Tt.An==A_int(ati),:);        
        Dens =100.*Tta.Aws./Tta.nwsr; Dens(isinf(Dens))=nan;
        R = Tta.Aws./Tta.At; 
        DataRWS(:,ati,t)= [nanmean(R) CI(R,prob)];
        Data(:,ati,t)=[nanmean(Dens) CI(Dens,prob)] ; 
    end
    
end


figure;
subplot(1,3,1)
hold on; grid on
x = A_int; Data1 = DataRWS;
xr = [x fliplr(x)];
for t=1:3 % Plot Data
    yr = [Data1(1,:,t)+Data1(2,:,t), fliplr(Data1(1,:,t)+Data1(3,:,t))];
    f=fill(xr,yr,cmap(t,:));
    set(f,'facealpha',.2,'LineStyle','none')
    p(t) = plot(x,Data1(1,:,t), 'Color',cmap(t,:), 'LineWidth', 2);
end
title('a) Ratio empty space to slime mold')
xlabel('Normalized area [-]'); ylabel('Empty space over slime mold areas [-]'); legend([p(1) p(2) p(3)],LegendTreatments,'Location', 'Best')

Tta = T(T.An==A_int(end),:);  

shortlegends = {'Neu','Nut','Adv'};
subplot(1,3,2);
boxplot(Tta.nwsr,Tta.Treatment,'orientation','vertical','label',shortlegends,'color',cmap,'OutlierSize',1); 
ylabel('Number of enclosed regions [-]'); grid on
title('b) Number of empty regions')

shortlegends = {'Neu','Nut','Adv'};
subplot(1,3,3);
boxplot(Tta.MeanEnclosed,Tta.Treatment,'orientation','vertical','label',shortlegends,'color',cmap,'OutlierSize',1); 
ylabel('Size of empty region [mm^2]'); grid on
title('c) Mean size of empty region')

% plotTreatmentsTime_2(A_int,DataRWS,Data,{'a) Ratio empty space to slime mold','b) Mean size of empty region',},'Normalized area [-]',...
%     {'Empty space over slime mold areas [-]','Size of empty region [mm^2]'},LegendTreatments)
% subplot(1,2,2); xlim([2 4])
% saveas(gcf,'BF_RatioEmptyArea.png')
% print('BF_RatioEmptyArea','-depsc');
% saveas(gcf,'BF_RatioEmptyArea.fig')

% plotTreatmentsTime(A_int,Data,'Mean size of enclosed region [cm^2]',...
%     'Normalized Area [-]','Mean size of enclosed region [cm^2]',Treatments)
% xlim([2 4])
% saveas(gcf,'BF_AreaEnclosedRegions.png')
% print('BF_AreaEnclosedRegions','-depsc');
% saveas(gcf,'BF_AreaEnclosedRegions.fig')
% 
% plotTreatmentsTime(A_int,DataRWS,'Ratio Empty Area to SM Area',...
%     'Normalized Area [-]','Area Ratio [-]',Treatments)
% %xlim([2 4])
% saveas(gcf,'BF_RatioEmptyArea.png')
% print('BF_RatioEmptyArea','-depsc');
% saveas(gcf,'BF_RatioEmptyArea.fig')

% 4) Percentage of Clusters

Data = zeros(3,numel(A_int),numel(Treatments));

for t=1:3 % Gather Data
    Tt = T(strcmp(T.Treatment,Treatments{t}),:);    
    for ati=1:numel(A_int)
        Tta = Tt(Tt.An==A_int(ati),:);
        pC =Tta.pC;
        Data(:,ati,t)=[mean(pC) CI(pC,prob)] ; 
    end    
end

figure;
subplot(1,2,1)
plotTreatmentsTime(A_int,Data,'Percentage of cell corresponding to Clusters',...
    'Normalized Area [-]','Cluster Area [%]',LegendTreatments)
subplot(1,2,2)
Tt = T(T.An==4,:);
boxplot(Tt.pC,Tt.Treatment,'orientation','vertical','label',LegendTreatments,'color',cmap); 
ylabel('Percentage of Clusters [%] at An=4'); grid on; %ylim([0 0.5])

saveas(gcf,'BF_pcntClusters.png')
% saveas(gcf,'BF_pcntClusters.fig')

%% 5) Cluster vs vein area

% DataC = zeros(3,numel(A_int),numel(Treatments));
% DataV = zeros(3,numel(A_int),numel(Treatments));
% 
% for t=1:3 % Gather Data
%     Tt = T(strcmp(T.Treatment,Treatments{t}),:);
%     
%     for ati=1:numel(A_int)
%         Tta = Tt(Tt.An==A_int(ati),:);
%         pC = Tta.An.*Tta.pC./100;
%         DataC(:,ati,t)=([mean(pC) CI(pC,p)]);
%         pC = Tta.An.*(100-Tta.pC)./100;
%         DataV(:,ati,t)=([mean(pC) CI(pC,p)]) ;
%     end
%     
% end
% 
% figure(5)
% hold on
% colorsCode = jet(3);
% nanvector = nan(1,numel(A_int));
% for t=1:3 % Plot Data
%     subplot(1,2,1)
%     hold on
%     errorbar(A_int,DataC(1,:,t),DataC(2,:,t),DataC(3,:,t),nanvector,nanvector)
%     subplot(1,2,2)
%     hold on
%     errorbar(A_int,DataV(1,:,t),DataV(2,:,t),DataV(3,:,t),nanvector,nanvector)
% end
% 
% subplot(1,2,1)
% ylabel('Clusters Area [-]')
% xlabel('Normalized Area [-]')
% 
% subplot(1,2,2)
% ylabel('Veins Area [-]')
% xlabel('Normalized Area [-]')
% legend(Treatments)
% 
% sgtitle('Area distribution Veins/Custers')
% saveas(gcf,'BF_VeinsVsClusters.png')

%% Mean W-L Distributions

Tta = T(T.An==A_int(end),:);   % Only at An=4

W_all = Tta(:,10:29);
L_all = Tta(:,30:49);

nBins=20;
DataW = zeros(3,nBins,numel(Treatments));
DataL = zeros(3,nBins,numel(Treatments));

for t=1:numel(Treatments)
    Tt = T(strcmp(Tta.Treatment,Treatments{t}),:);
    W = table2array(W_all(strcmp(Tta.Treatment,Treatments{t}),:));
    L = table2array(L_all(strcmp(Tta.Treatment,Treatments{t}),:));
     
    DataW(:,:,t)=100.*[mean(W); CI(W,prob,2)] ; 
    DataL(:,:,t)=100.*[mean(L); CI(L,prob,2)] ;
end

Wbins=linspace(0,2.5,nBins+1); % Width Bins -> Increase!
Lbins=linspace(0,6,nBins+1); % Length Bins
meanW=(Wbins(2:end)+Wbins(1:end-1))./2;
meanL=(Lbins(2:end)+Lbins(1:end-1))./2;

figure; % Width-Length Distributions

for t=1:3 % Plot Data
    subplot(2,3,t)
    histogram('BinCounts', DataW(1,:,t), 'BinEdges', Wbins,'FaceAlpha',0.4,'FaceColor',cmap(t,:),'EdgeColor','none')
    
    hold on
    er = errorbar(meanW,DataW(1,:,t),DataW(1,:,t)-DataW(2,:,t));    
    er.Color = cmap(t,:); er.LineStyle = 'none';  
    ylabel('Veins frequency [%]'); xlabel('Vein Width [mm]')
    yticks(0:2.5:22.5)
    ylim([0 22.5]); xlim([0 2.5]); title(Treatments{t}); grid on

    subplot(2,3,t+3)
    histogram('BinCounts', DataL(1,:,t), 'BinEdges', Lbins,'FaceAlpha',0.4,'FaceColor',cmap(t,:),'EdgeColor','none')
    
    hold on
    er = errorbar(meanL,DataL(1,:,t),DataL(1,:,t)-DataL(2,:,t)); er.Color = cmap(t,:);                            
    er.LineStyle = 'none';  
    
    ylabel('Veins frequency [%]'); xlabel('Vein Length [%]')
    yticks(0:2.5:17.5)
    ylim([0 17.5]); xlim([0 6]); grid on ;% title(Treatments{t})
end

%sgtitle('Distribution of vein length and width by Treatment')


a=1;


%% 6) W-L distribution: Difference among treatments

% comp=[2 1;3 1;2 3];
% comp_str = {'Glucose-Control','NaCl - Control','Glucose - NaCl'};
% BTitles={'Diff between substrate - Before Fusion'};
% 
% At = 4;
% %pxpmm = 1100/90; % this many pixels corresponds to 1mm.
% 
% nBins=20;
% DataW = zeros(3,nBins,numel(Treatments));
% Data = zeros(3,nBins,numel(Treatments));
% 
% SampleExp = zeros(2,nBins,numel(Treatments)); % First Row W - Second Row L
% 
% Wbins=linspace(0,2.5,nBins+1); % Width Bins -> Increase!
% Lbins=linspace(0,6,nBins+1); % Length Bins
% meanW=(Wbins(2:end)+Wbins(1:end-1))./2;
% meanL=(Lbins(2:end)+Lbins(1:end-1))./2;
% 
% for t=1:numel(Treatments)
%     
%     Tt = T(strcmp(T.Treatment,Treatments{t}),:);
%     
%     Tta = Tt(Tt.An>=2,:);
%     W = Tta.Wbins;
%     L = Tta.Lbins;
%     DataW(:,:,t)=[mean(W); CI(W,prob,2)] ;
%     Data(:,:,t)=[mean(L); CI(L,prob,2)] ;
%     
%     Tta = Tt(Tt.An==At,:);
%     
%     nSample =3;
%     SampleExp(1,:,t) = Tta.Wbins(nSample,:);
%     SampleExp(2,:,t) = Tta.Lbins(nSample,:);    
%     
%     figure(8)
%     subplot(2,3,t)
%     hold on
%     histogram('BinCounts', SampleExp(1,:,t), 'BinEdges', Wbins,'FaceAlpha',0.55,'EdgeColor','none','Normalization','pdf') % 
%     pd = makedist('Loglogistic','mu',Tta.Wpar(2,1),'sigma',Tta.Wpar(2,2));
%     plot(Wbins,pdf(pd,Wbins))
%     title(Treatments{t})
%     
%     subplot(2,3,t+3)
%     hold on
%     histogram('BinCounts', SampleExp(2,:,t), 'BinEdges', Lbins,'FaceAlpha',0.55,'EdgeColor','none','Normalization','pdf') % 
%     pd = makedist('Loglogistic','mu',Tta.Lpar(2,1),'sigma',Tta.Lpar(2,2));
%     plot(Lbins,pdf(pd,Lbins))
%     
% end
% 
% sgtitle('Sample of veins Width and Length distributions, fitted to loglogistic')
% saveas(gcf,'BF_SampleWLdist.png')
% 
% figure(6) % Width Distribution
% 
% for t=1:3 % Plot Data
%     subplot(1,3,t)
%     histogram('BinCounts', DataW(1,:,t), 'BinEdges', Wbins,'FaceAlpha',0.75,'EdgeColor','none','Normalization','probability')
%     
%     hold on
%     er = errorbar(meanW,DataW(1,:,t),DataW(2,:,t),DataW(3,:,t));    
%     er.Color = [0 0 1];                            
%     er.LineStyle = 'none';  
% 
%     ylabel('Percentage of Veins [%]')
%     xlabel('Vein Width [mm]')
%     ylim([0 0.4])
%     title(Treatments{t})
% 
% end
% sgtitle('Distribution of vein width by Treatment')
% saveas(gcf,'BF_Wdist.png')
% 
% figure; % Length Distribution
% for t=1:3 % Plot Data
%     subplot(1,3,t)
%     histogram('BinCounts', Data(1,:,t), 'BinEdges', Lbins,'FaceAlpha',0.5,'EdgeColor','none','Normalization','probability')
%     
%     hold on
%     er = errorbar(meanL,Data(1,:,t),Data(2,:,t),Data(3,:,t));    
%     er.Color = [0 0 1];                            
%     er.LineStyle = 'none';  
%     
%     ylabel('Percentage of Veins [%]')
%     xlabel('Vein Length [%]')
%     ylim([0 0.25])
%     title(Treatments{t})
% 
% end
% sgtitle('Distribution of vein length by Treatment')
% saveas(gcf,'BF_Ldist.png')

%% 7) W-L distribution: Distribution of fitted parameters

Tta = T(T.An==2,:);
cat = Tta.Treatment; 
WData = [Tta.Wpar_1 Tta.Wpar_2];
LData = [Tta.Lpar_1 Tta.Lpar_2];

figure;
h = scatterhist(WData(:,1),WData(:,2),'Group',cat,'Marker','+o*','Color',cmap);
grid on
xlabel('mu [-]')
ylabel('theta [-]')
hold on;
boxplot(h(2),WData(:,1),cat,'orientation','horizontal','label',LegendTreatments,'color',cmap);
grid(h(2),'on')
boxplot(h(3),WData(:,2),cat,'orientation','horizontal','label', LegendTreatments,'color',cmap);
grid(h(3),'on')
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
axis(h(1),'auto');  % Sync axes
hold off;
title('Width distribution: Fitted parameters')
% saveas(gcf,'BF_Wfit.png')
% print('BF_Wfit','-depsc');
% 
% figure;
% h = scatterhist(LData(:,1),LData(:,2),'Group',cat,'Marker','+o*','Color',cmap);
% grid on
% xlabel('mu [-]')
% ylabel('theta [-]')
% hold on;
% boxplot(h(2),LData(:,1),cat,'orientation','horizontal','label',LegendTreatments,'color',cmap);
% grid(h(2),'on')
% boxplot(h(3),LData(:,2),cat,'orientation','horizontal','label', LegendTreatments,'color',cmap);
% grid(h(3),'on')
% set(h(2:3),'XTickLabel','');
% view(h(3),[270,90]);  % Rotate the Y plot
% axis(h(1),'auto');  % Sync axes
% hold off;
% title('Length distribution: Fitted parameters')
% saveas(gcf,'BF_Lfit.png')
% saveas(gcf,'BF_Lfit.fig')
% print('BF_Lfit','-depsc');


figure; 
shortlegends = {'Neu','Nut','Adv'};
subplot(1,4,1);
boxplot(T.Wpar_1,T.Treatment,'orientation','vertical','label',shortlegends,'color',cmap,'OutlierSize',1); 
ylabel('\mu '); grid on
title('Vein width distribution')
subplot(1,4,2);
boxplot(T.Wpar_2,T.Treatment,'orientation','vertical','label',shortlegends,'color',cmap,'OutlierSize',1); 
ylabel('\theta'); grid on; ylim([0.2 0.9])

subplot(1,4,3);
boxplot(T.Lpar_1,T.Treatment,'orientation','vertical','label',shortlegends,'color',cmap,'OutlierSize',1); 
ylabel('\mu '); grid on
title('Vein length distribution')
subplot(1,4,4);
boxplot(T.Lpar_2,T.Treatment,'orientation','vertical','label',shortlegends,'color',cmap,'OutlierSize',1); 
ylabel('\theta'); grid on; ylim([0.2 0.9])



% 8) Total Length of skeleton

Data = zeros(3,numel(A_int),numel(Treatments));
DataV = zeros(3,numel(A_int),numel(Treatments));
%A0 = (13/20)^2*pi(); %cm2 area of initial cell
for t=1:3 % Gather Data
    Tt = T(strcmp(T.Treatment,Treatments{t}),:);    
    for ati=1:numel(A_int)
        Tta = Tt(Tt.An==A_int(ati),:);
        L = Tta.Le_p_2; %(:,2);
        V = Tta.meanVein; V(isinf(V))=[];
        Data(:,ati,t)=[mean(L) CI(L,prob)] ;
        DataV(:,ati,t)=[nanmean(V) CI(V,prob)] ;
    end
end

% plotTreatmentsTime_2(A_int,Data,DataV,{'a) Network length','b) Vein width'},'Normalized Area [-]',...
%     {'Total network length [cm]','Average vein width [mm]'},Treatments)


figure;
subplot(1,3,1)
hold on; grid on
x = A_int; Data1 = Data;
xr = [x fliplr(x)];
for t=1:3 % Plot Data
    yr = [Data1(1,:,t)+Data1(2,:,t), fliplr(Data1(1,:,t)+Data1(3,:,t))];
    f=fill(xr,yr,cmap(t,:));
    set(f,'facealpha',.2,'LineStyle','none')
    p(t) = plot(x,Data1(1,:,t), 'Color',cmap(t,:), 'LineWidth', 2);
end
title('a) Network length')
xlabel('Normalized area [-]'); ylabel('Total network length [cm]'); legend([p(1) p(2) p(3)],LegendTreatments,'Location', 'Best')


subplot(1,3,2)
hold on; grid on
x = A_int; Data1 = DataV;
xr = [x fliplr(x)];
for t=1:3 % Plot Data
    yr = [Data1(1,:,t)+Data1(2,:,t), fliplr(Data1(1,:,t)+Data1(3,:,t))];
    f=fill(xr,yr,cmap(t,:));
    set(f,'facealpha',.2,'LineStyle','none')
    p(t) = plot(x,Data1(1,:,t), 'Color',cmap(t,:), 'LineWidth', 2);
end
title('b) Vein width')
xlabel('Normalized area [-]'); ylabel('Average vein width [mm]'); legend([p(1) p(2) p(3)],LegendTreatments,'Location', 'Best')


Tta = T(T.An==A_int(end),:);  

shortlegends = {'Neu','Nut','Adv'};
subplot(1,3,3);
boxplot(Tta.meanVein,Tta.Treatment,'orientation','vertical','label',shortlegends,'color',cmap,'OutlierSize',1); 
ylabel('Average vein width [mm]'); grid on
title('c) Vein Width at An=4')
ylim([0.4 1.4])


% subplot(1,2,2)
% xlim([2 4])
% plotTreatmentsTime(A_int,Data,'Total network length',...
%     'Normalized Area [-]','Skeleton Length [cm]',Treatments)
% saveas(gcf,'BF_NetworkLength.png')
% saveas(gcf,'BF_NetworkLength.fig')
% 
% plotTreatmentsTime(A_int,DataV,'Average Vein thickness',...
%     'Normalized Area [-]','Average Veins Thickness [mm]',Treatments)
% xlim([2 4])

% saveas(gcf,'BF_Length_Width.png')
% saveas(gcf,'BF_Length_Width.fig')
% print('BF_Length_Width','-depsc');


writetable(T,fullfile(DataFolder,'Data_BF.xlsx'),'Sheet',1,'Range','A1')

% for t=1:3 % Combinations
% 
%     % W Bins
%     Tta = Tt(Tt.An==At),:);
%     pC =Tta.pC;
%         
%      sum(T.Lbins,2);
%     figure(5+ba) %% Difference Before fusion
%     subplot(2,3,t) % width
%     hold on
%     yval=100*(HistData{ba,comp(t,1),1}-HistData{ba,comp(t,2),1});    
%     positiveIndexes = yval >= 0; % Find where data is positive or negative.
%     negativeIndexes = yval <  0;
%     bar(bins{1}(positiveIndexes), yval(positiveIndexes), 'b', 'BarWidth', 1)
%     bar(bins{1}(negativeIndexes), yval(negativeIndexes), 'r', 'BarWidth', 1)
%     xlabel('Width [px]')
%     ylabel('Difference between treatments [%]')
%     title(comp_str{t})
% 
%     subplot(2,3,3+t) % Length
%     hold on
%     yval=100*(HistData{ba,comp(t,1),2}-HistData{ba,comp(t,2),2});    
%     positiveIndexes = yval >= 0; % Find where data is positive or negative.
%     negativeIndexes = yval <  0;
%     bar(bins{2}(positiveIndexes), yval(positiveIndexes), 'b', 'BarWidth', 1)
%     bar(bins{2}(negativeIndexes), yval(negativeIndexes), 'r', 'BarWidth', 1)
%     xlabel('Length [px]')
%     ylabel('Difference between treatments [%]')   
% end
% sgtitle(BTitles{ba});




% 
% save(fullfile(DataFolder,Folders{f},'PreFusion_indexes.mat'),'Data_f','-v7.3')
% 
%     save(fullfile(DataFolder,Folders{f},'AfterFusion_indexes.mat'),'Data_f','-v7.3');      
%     
%     %Save Fusion Regions
%     save(fullfile(DataFolder,Folders{f},'FusionRegions.mat'),'FusionRegionsExperiments','-v7.3');    
%     
%     
%      if (isfile(fullfile(DataFolder,Folders{f},srcFiles_filename))==2) % If exists load it, otherwise compute it.
%         srcFiles = load(fullfile(DataFolder,Folders{f},srcFiles_filename));



%%% AUXILIAR FUNCTIONS %%%


function [names]=GetSubfolders(path)
        d = dir(path);
        isub = [d(:).isdir]; %# returns logical vector
        names = {d(isub).name}';
        names(ismember(names,{'.','..'})) = [];
end 

function [I]=CI(x,p,dim)
    
if isvector(x)    
     I = nanstd(x)/sqrt(sum(~isnan(x))) * tinv(abs([0,1]-(1-p)/2),sum(~isnan(x))-1) ;%+ nanmean(x);    
    pp=1;
    
else
    if dim == 2 % columns
        nC = size(x,2);
        I = zeros(2,nC);   
            for c=1:nC
                I(:,c) = (nanstd(x(:,c))/sqrt(sum(~isnan(x(:,c)))) * tinv(abs([0,1]-(1-p)/2),sum(~isnan(x(:,c)))-1) + nanmean(x(:,c)))';    
            end    
    elseif dim==1 % rows
            nR = size(x,dim);
            I = zeros(nR,2);   
            for r=1:nR
                I(:,r) = nanstd(x(r,:))/sqrt(nansum(x(r,:))) * tinv(abs([0,1]-(1-p)/2),nansum(x(r,:))-1) + nanmean(x(r,:));    
            end  

    end
end


end

function[]=plotTreatmentsTime(x,Data,tit,xl,yl,leg)

cmap = [27,158,119; 217,95,2; 117,112,179]./255;

%figure; 
hold on

xr = [x fliplr(x)];
for t=1:3 % Plot Data
    yr = [Data(1,:,t)+Data(2,:,t), fliplr(Data(1,:,t)+Data(3,:,t))];
    f=fill(xr,yr,cmap(t,:));
    set(f,'facealpha',.2,'LineStyle','none')
    p(t) = plot(x,Data(1,:,t), 'Color',cmap(t,:), 'LineWidth', 2);
end

xlabel(xl); ylabel(yl); title(tit); 
legend([p(1) p(2) p(3)],leg,'Location', 'Best')
%saveas(gcf,[savepath '.png'])

end

function[]=plotTreatmentsTime_2(x,Data1,Data2,tit,xl,yls,leg)

cmap = [27,158,119; 217,95,2; 117,112,179]./255;

figure;
subplot(1,2,1)
hold on
grid on
xr = [x fliplr(x)];
for t=1:3 % Plot Data
    yr = [Data1(1,:,t)+Data1(2,:,t), fliplr(Data1(1,:,t)+Data1(3,:,t))];
    f=fill(xr,yr,cmap(t,:));
    set(f,'facealpha',.2,'LineStyle','none')
    p(t) = plot(x,Data1(1,:,t), 'Color',cmap(t,:), 'LineWidth', 2);
end
title(tit{1})
xlabel(xl); ylabel(yls{1}); legend([p(1) p(2) p(3)],leg,'Location', 'Best')

subplot(1,2,2)
hold on
grid on
xr = [x fliplr(x)];
for t=1:3 % Plot Data
    yr = [Data2(1,:,t)+Data2(2,:,t), fliplr(Data2(1,:,t)+Data2(3,:,t))];
    f=fill(xr,yr,cmap(t,:));
    set(f,'facealpha',.2,'LineStyle','none')
    p(t) = plot(x,Data2(1,:,t), 'Color',cmap(t,:), 'LineWidth', 2);
end
title(tit{2})
xlabel(xl); ylabel(yls{2}); % legend([p(1) p(2) p(3)],leg,'Location', 'Best')
%sgtitle(tit);

%saveas(gcf,[savepath '.png'])

end

