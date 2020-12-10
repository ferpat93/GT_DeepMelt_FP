close all

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.065 0.065], [0.2 0.05], [0.075 0.04]); % subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
if ~make_it_tight,  clear subplot;  end

%% ~ AFTER - FUSION ~ %%

% Constants
cmap = [27,158,119; 217,95,2; 117,112,179]./255; % Colorblind safe colormap
prob = 0.95; % Confidence interval
TW = 180; % Time window (before/after fusion)
dt = 15; % Time Interval
timesAF = (0:dt:TW)'; Ao_mm = (13/20)^2*pi(); %mm2
%A_int = 1:0.5:4; % Area Intervals for part 1 (relative to initial area)
DishTreatments = {'Control','Control','Glucose','Glucose','NaCl','NaCl',...
    'Control','Control','Glucose','Glucose','NaCl','NaCl',...
    'Glucose','Glucose','Glucose','Glucose',...
    'NaCl','NaCl','NaCl','NaCl'};

Treatments = {'Control','Glucose','NaCl'};
shortlegends = {'Neu','Nut','Adv'};
LegendTreatments = {'Neutral','Nutritive','Adverse'};

DataFolder = 'D:\Slime_Mold_Network\Data\Fusion';
DataFolder = '/Users/lfp3/Dropbox (GaTech)/GT/Spring-20/SMN/Data';
ResultsFolder = 'D:\Slime_Mold_Network\Results\Fusion';

LedgerName = 'Ledger.xlsx';
T = readtable(fullfile(DataFolder,LedgerName));

% Gather Data

Folders=GetSubfolders(DataFolder);
reload = false;

if or(~isfile(fullfile(DataFolder,'Data_SMN_AF_3h.mat')),reload)
    
    TableVariables = {'nDish','tfusion','t','Ao','At','pC','Ae','Aws','nwsr','Wbins','Lbins','Wpar','Lpar','Inds','folder'};
    varTypes = [repelem({'double'},numel(TableVariables)-1)  {'string'}];
    sz = [0 numel(TableVariables)];
    TG = table('Size',sz,'VariableTypes',varTypes,'VariableNames',TableVariables); % General Indexes table

    TableVariables = { 'Ao' 'Awso' 'Lo' 'At' 'Aws' 'Ac' 'Lt' 'Lv' 'Dt' 'Dv','maxD','t',...
        'nNodes_tot','nEdges_tot','nNodes_in','nEdges_in','n_Edges_cross','length_path','drag_path','mfn','mfw',...
        'Treat','Folder'}; % Fusion Region Indexes table
    varTypes = [repelem({'double'},numel(TableVariables)-2)  {'string'} {'string'}];
    TF = table('Size',[0 numel(TableVariables)],'VariableTypes',varTypes,'VariableNames',TableVariables);
    
    for f=1:numel(Folders)

        if isfile(fullfile(DataFolder,Folders{f},'AfterFusion_indexes.mat'))
            Data_fi = extractfield(load(fullfile(DataFolder,Folders{f},'AfterFusion_indexes.mat')),'Data_f');
            TG = vertcat(TG,Data_fi{1});
        end

        if isfile(fullfile(DataFolder,Folders{f},'FusionRegions.mat'))
            Regions = extractfield(load(fullfile(DataFolder,Folders{f},'FusionRegions.mat')),'FusionRegionsExperiments');
            Regions = Regions{1};

            for i=1:size(Regions,1)
                nDish = i;
                treat = DishTreatments{nDish};          

                if ~isempty(Regions{nDish,1})
                    ti = getFRindexes(Regions{nDish,1},treat,Folders{f},timesAF);    
                    [tg, ppmm] = getGraphIndexes(Regions(nDish,:),timesAF);
                    KA = (ppmm)^2;
                    ti.Ao = ti.Ao.*KA; ti.Awso = ti.Awso.*KA; ti.At = ti.At.*KA; ti.Aws = ti.Aws.*KA;
                    ti.Dt = ti.Dt.*ppmm; ti.Dv = ti.Dv.*ppmm; 
                    ti.Lo = ti.Lo.*ppmm; ti.Lt = ti.Lt.*ppmm; ti.Lv = ti.Lv.*ppmm; ti.maxD = ti.maxD.*ppmm;
                    TF = [TF ; horzcat(tg,ti)];
                    
                end
            end
        end

    end

    TG = addvars(TG,DishTreatments(TG.nDish)','NewVariableNames','Treatment');
    TG(TG.t==inf,:)=[]; % Erase null entries
    TG(isnan(TG.Ao),:)=[]; % Erase null entries

    save(fullfile(DataFolder,'Data_SMN_AF_3h.mat'),'TG','TF','-v7.3');
    writetable(TG,fullfile(DataFolder,'Data_AF.xlsx'),'Sheet',1,'Range','A1')
    writetable(TF,fullfile(DataFolder,'Data_FusionRegion.xlsx'),'Sheet',1,'Range','A1')
else
    load(fullfile(DataFolder,'Data_SMN_AF_3h.mat'),'TG','TF');
end


TF.nEdges_in_cross = TF.nEdges_in+TF.n_Edges_cross;
TF.ID = repelem((1:81)',13);
TF.AN = TF.At./TF.Ao;
TF.LN = TF.Lt./TF.Lo;
TF.PAN = (TF.At + TF.Aws)./(TF.Ao + TF.Awso);
TF.ES = TF.Aws./TF.At;
NET = TF.nEdges_in+TF.n_Edges_cross; NEI = TF.nEdges_in; NN = TF.nNodes_in;
TF.AEI = (NEI-NN+1)./(2*NEI-5); 
TF.AET = (NET-NN+1)./(2*NET-5); 
TF.AT = (TF.nEdges_tot-TF.nNodes_tot+1)./(2*TF.nEdges_tot-5);
TF.mfmw = TF.mfw./TF.mfn;

% save(fullfile(DataFolder,'Data_SMN_AF_3h.mat'),'TG','TF','-v7.3');
% writetable(TG,fullfile(DataFolder,'Data_AF.xlsx'),'Sheet',1,'Range','A1')
writetable(TF,fullfile(DataFolder,'Data_FusionRegion.xlsx'),'Sheet',1,'Range','A1')
a=1;

    

%% -- START PLOTTING -- %%
 

% RENTIAN SCALING
%gscatter(log(TF.nNodes_in),log(TF.n_Edges_cross),TF.Treat)


%% Nodes vs edges

% figure; 
% %subplot(1,2,1);
% hold on;
% 
% NET = TF.nEdges_in+TF.n_Edges_cross; NEI = TF.nEdges_in; NN = TF.nNodes_in;
% 
% m=NN\NEI; rsq=1-sum((NEI-NN*m).^2)/sum((NEI-mean(NEI)).^2);
% line([0 max(NN)],m.*[0 max(NN)],'Color',[0.5 0 0],'LineStyle','--')
% txt1 = ['EIR = ' num2str(m,3) 'NR - R^2 = ' num2str(rsq,2)]; t2=text(0.4, 0.25, txt1,'Units','normalized');
% set(t2,'Rotation',30)
% m=NN\NET; rsq=1-sum((NET-NN*m).^2)/sum((NET-mean(NET)).^2);
% line([0 max(NN)],m.*[0 max(NN)],'Color',[0.5 0 0],'LineStyle','-.')
% txt1 = ['EAR = ' num2str(m,3) 'NR - R^2 = ' num2str(rsq,2)]; t1=text(0.4, 0.4, txt1,'Units','normalized');
% set(t1,'Rotation',37)
% xlim([0 200]); ylim([0 350])
% xlabel('Number of Nodes - FR'); ylabel('Number of Edges - FR')
% legend({'Edges Inside FR','Edges Across FR'},'Location','northwest');
% 
% ax1 = gca; % current axes
% ax1.XColor = [0.5 0 0]; ax1.YColor = [0.5 0 0];
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');
% 
% m=TF.nNodes_tot\TF.nEdges_tot; rsq=1-sum((TF.nEdges_tot-TF.nNodes_tot*m).^2)/sum((TF.nEdges_tot-mean(TF.nEdges_tot)).^2);
% line([0 max(TF.nNodes_tot)],m.*[0 max(TF.nNodes_tot)])
% txt1 = ['ECN = ' num2str(m,3) 'NN - R^2 = ' num2str(rsq,2)]; t3=text(0.7, 0.575, txt1,'Units','normalized');
% set(t3,'Rotation',33.5)
% xlim([0 3000]); ylim([0 5250])
% xlabel('Number of Nodes - CN'); ylabel('Number of Edges - CN')
% legend('Complete Network');

%saveas(gcf,'AF_EdgesVsNodes.png');  saveas(gcf,'AF_EdgesVsNodes.fig');

%% Boxplot
% alpha = [TF.AEI TF.AET TF.AT];
% data = {alpha(strcmp(TF.Treat,Treatments{1}),:),alpha(strcmp(TF.Treat,Treatments{2}),:),alpha(strcmp(TF.Treat,Treatments{3}),:)};
% 
% subplot(1,2,2);
% boxplotGroup(data, 'PrimaryLabels', {'C','G','N'}, ...
%   'SecondaryLabels',{'Inside FR', 'Across FR','Complete Network'}, 'GroupLabelType', 'Horizontal')
% grid on
% title('Connectivity of fusion region vs. complete network')
% ylim([0 0.5]); ylabel('Alpha')

% saveas(gcf,'AF_GraphsConnectivity.png');  saveas(gcf,'AF_GraphsConnectivity.fig');
% print('AF_GraphsConnectivity','-depsc');

%% Connecting edges and "bandwidth"
DataN = zeros(3,numel(timesAF),numel(Treatments));
DataW = zeros(3,numel(timesAF),numel(Treatments));
DataT = zeros(3,numel(timesAF),numel(Treatments));
prob = 0.95;
for t=1:3 % Gather Data
    Tt = TF(strcmp(TF.Treat,Treatments{t}),:);
    %Tt(isinf(Tt.length_path),:)=[];
    for ati=1:numel(timesAF)
        Tta = Tt(Tt.t==timesAF(ati),:);
        N = Tta.mfn;
        W = Tta.mfw;
        T = W./N;
        DataN(:,ati,t)=[nanmean(N) CI(N,prob)] ; 
        DataW(:,ati,t)=[nanmean(W) CI(W,prob)] ; 
        DataT(:,ati,t)=[nanmean(T) CI(T,prob)] ; 
    end   
end

%plotBigraphTime(timesAF/60,DataN,DataW,'Path change between initial cells','Time from fusion [h]',...
%    {'Cells Connectivity after Fusion','BandWidth [mm]'},Treatments)
Data4{1}=DataN; Data4{2}=DataT; 
% plotBigraphTime(timesAF/60,DataN,DataT,'Path change between initial cells','Time from fusion [h]',...
%     {'Cells Connectivity after Fusion','Mean Vein width [mm]'},Treatments)
% saveas(gcf,'AF_ConnectingVeins.png');  saveas(gcf,'AF_ConnectingVeins.fig');

% TT = TF(TF.t==timesAF(end),:);
% nv = TT.mfn;
% mw = TT.mfw./TT.mfn;
% figure;
% subplot(1,2,1)
% boxplot(nv,TT.Treat,'orientation','vertical','label',Treatments,'color',cmap);grid on
% xlabel('Treatments'); ylabel('Number of connecting veins')
% subplot(1,2,2)
% boxplot(mw,TT.Treat,'orientation','vertical','label',Treatments,'color',cmap);grid on
% xlabel('Treatments'); ylabel('Average connecting vein width')
% sgtitle('Connecting Veins 3h After Fusion')
% saveas(gcf,'ConnVeins_After3h.png')

%% Evolution of L and D

ind = find(TF.t==0);
DataD = zeros(3,numel(timesAF),numel(Treatments));
DataL = zeros(3,numel(timesAF),numel(Treatments));
prob = 0.95;
for t=1:3 % Gather Data
    Tt = TF(strcmp(TF.Treat,Treatments{t}),:);
    Tt(isinf(Tt.length_path),:)=[];
    for ati=1:numel(timesAF)
        Tta = Tt(Tt.t==timesAF(ati),:);
        D = Tta.length_path;
        L = Tta.drag_path;
        DataD(:,ati,t)=[nanmean(D) CI(D,prob)] ; 
        DataL(:,ati,t)=[nanmean(L) CI(L,prob)] ; 
    end
    
end

Data4{3}=DataD; Data4{4}=DataL; 

% plot4Time(timesAF/60,Data4,{'a) Connecting veins between cells','b) Connecting veins width','c) Path length between cells','d) Drag between cells'},...
%     'Time from fusion [h]',...
%     {'Number of connecting veins','Average vein width [mm]','Path length [mm]','Path drag [mm^{-3}]'},...
%     LegendTreatments)

% 
% plotBigraphTime(timesAF/60,DataN,DataT,'Path change between initial cells','Time from fusion [h]',...
%     {'Cells Connectivity after Fusion','Mean Vein width [mm]'},Treatments)
% saveas(gcf,'AF_ConnectingVeins.png');  saveas(gcf,'AF_ConnectingVeins.fig');
% 
% plotBigraphTime(timesAF/60,DataD,DataL,'Path change between initial cells','Time from fusion [h]',...
%     {'Path Length [-]','Path Drag [-]'},Treatments)
% saveas(gcf,'AF_PathBetweenCells.png');  saveas(gcf,'AF_PathBetweenCells.fig');

% T1 = TF(TF.t==timesAF(1),:); T4 = TF(TF.t==timesAF(end-4),:);
% one = T1.drag_path;
% four = T4.drag_path;
% 
% figure;
% boxplot(100.*(four-one)./one,T1.Treat)


% 
% figure;
% for i=1:numel(ind)
%     figure(111); hold on
%     plot(timesAF,TF.length_path(ind(i):ind(i)+12))
%     
%     figure(112); hold on
%     %plot(times,TF.drag_path(ind(i):ind(i)+12)./TF.drag_path(ind(i)))
%     plot(timesAF,TF.drag_path(ind(i):ind(i)+12))
%     
% end

%% 1) Fusion -> Time and Area

F = TG(TG.tfusion==0,:); % Table with just entries at fusion

cat = F.Treatment; 
Data = [F.t./60 F.At./F.Ao]; % [time Area(norm)]

% figure; 
% subplot(1,2,1);
% boxplot(F.t./60,cat,'orientation','vertical','label',LegendTreatments,'color',cmap); 
% ylabel('Time to fusion [h]'); grid on; title('a) Time to fusion') %ylim([0 0.5])
% subplot(1,2,2);
% boxplot(F.At./F.Ao,cat,'orientation','vertical','label',LegendTreatments,'color',cmap); 
% ylabel('Normalized cell area at fusion [-]'); grid on; title('b) Area at fusion') %ylim([0 0.5])

%ScatterHistogram(Data,cat,{'Time to fusion [h]','Normalized area at fusion [-]','Area and time to fusion'},0)
%saveas(gcf,'AF_ATtoFusion.png');  saveas(gcf,'AF_ATtoFusion.fig');

%% 2) Change of length and area over time - FusionRegion

DataA = zeros(3,numel(timesAF),numel(Treatments));
DataL = zeros(3,numel(timesAF),numel(Treatments));

for t=1:3 % Gather Data
    Tt = TF(strcmp(TF.Treat,Treatments{t}),:);
    
    for ati=1:numel(timesAF)
        Tta = Tt(Tt.t==timesAF(ati),:);
        A = Tta.At./Tta.Ao;
        L = Tta.Lt./Tta.Lo;
        DataA(:,ati,t)=[mean(A) CI(A,prob)] ; 
        DataL(:,ati,t)=[mean(L) CI(L,prob)] ; 
    end
    
end
 
Data6{1}=DataA; Data6{2}=DataL;

% plotBigraphTime(timesAF/60,DataA,DataL,'Change of area and network length over time','Time from fusion [h]',...
%     {'Normalized Area [-]','Normalized skeleton length [-]'},Treatments)

%saveas(gcf,'AF_AnLoverTime.png'); saveas(gcf,'AF_AnLoverTime.fig')

%% 3) Max and mean Dist

DataM = zeros(3,numel(timesAF),numel(Treatments)); % max
DataA = zeros(3,numel(timesAF),numel(Treatments)); % average

for t=1:3 % Gather Data
    Tt = TF(strcmp(TF.Treat,Treatments{t}),:);
    
    for ati=1:numel(timesAF)
        Tta = Tt(Tt.t==timesAF(ati),:);
        Max = Tta.maxD;
        Avg = Tta.Dt;
        DataM(:,ati,t)=[nanmean(Max) CI(Max,prob)] ; 
        DataA(:,ati,t)=[nanmean(Avg) CI(Avg,prob)] ; 
    end
    
end

Data6{3}=DataA; Data6{4}=DataM;

% plotBigraphTime(timesAF/60,DataA,DataM,'Width change over time','Time from fusion [h]',...
%     {'Average width [-]','Maximum width [-]'},Treatments)

%saveas(gcf,'AF_MeanMaxVeinWidth.png'); saveas(gcf,'AF_MeanMaxVeinWidth.fig');

%% 4) Solidity and enclosed area

DataS = zeros(3,numel(timesAF),numel(Treatments)); % Solidity
DataEA = zeros(3,numel(timesAF),numel(Treatments)); % Enclosed Area

for t=1:3 % Gather Data
    Tt = TF(strcmp(TF.Treat,Treatments{t}),:);    
    for ati=1:numel(timesAF)
        Tta = Tt(Tt.t==timesAF(ati),:);
        EA = (Tta.At + Tta.Aws)./(Tta.Ao + Tta.Awso);
        S = Tta.Aws./Tta.At;
        DataEA(:,ati,t)=[nanmean(EA) CI(EA,prob)] ; 
        DataS(:,ati,t)=[nanmean(S) CI(S,prob)] ; 
    end
end
 
Data6{5}=DataEA; Data6{6}=DataS;

plot3Time(timesAF/60,Data6([1 5 6]),{'Slime mold area','Print area','Empty area - slime mold ratio'},...
    'Time from fusion [h]',...
    {'Normalized slime mold area [-]','Normalized print area [-]','Empty area over slime mold [-]'},LegendTreatments)

plot3Time(timesAF/60,Data6([2 3 4]),{'Normalized network length','Average vein width','Maximum vein width'},...
    'Time from fusion [h]',...
    {'Normalized network length [-]','Average vein width [mm]','Maximum vein width [mm]'},LegendTreatments)


% plotBigraphTime(timesAF/60,DataA,DataL,'Change of area and network length over time','Time from fusion [h]',...
%     {'Normalized Area [-]','Normalized skeleton length [-]'},Treatments)
% plotBigraphTime(timesAF/60,DataA,DataM,'Width change over time','Time from fusion [h]',...
%     {'Average width [-]','Maximum width [-]'},Treatments)
% plotBigraphTime(timesAF/60,DataEA,DataS,'Slime Mold Print Evolution','Time from fusion [h]',...
%     {'Print Area [-]','Solidity [-]'},Treatments)

%saveas(gcf,'AF_EnclosedArea.png'); saveas(gcf,'AF_EnclosedArea.fig')

%% 5) Change of width and Length distributions

% time = 180; % time before/after fusion
% nBins=20;
% Wbins=linspace(0,2.5,nBins+1); % Width Bins 
% Lbins=linspace(0,6,nBins+1); % Length Bins
% nanvector = nan(1,nBins);
% 
% sumdiff = [];
% cat ={};
% for t=1:3 % Gather Data
%     
%     Tt = TG(strcmp(TG.Treatment,Treatments{t}),:); 
%     TtA = Tt(Tt.tfusion==time,:); % After
%     TtB = Tt(Tt.tfusion==-1.*time,:); % Before
%     
%     DataB = zeros(size(TtB,1),nBins,2); % Front: Width / Back: Length
%     DataA = zeros(size(TtA,1),nBins,2); % Front: Width / Back: Length
% 
%     DataA (:,:,1) = TtA.Wbins;
%     DataA (:,:,2) = TtA.Lbins;
%     
%     DataB (:,:,1) = TtB.Wbins;
%     DataB (:,:,2) = TtB.Lbins;
%     
%     [~,ia] = intersect([TtB.folder TtB.Treatment TtB.nDish],[TtA.folder TtA.Treatment TtA.nDish],'rows');
%     ChangeW = TtA.Wbins-TtB.Wbins(ia); 
%     ChangeL = TtA.Lbins-TtB.Lbins(ia);
% 
%     DataW=[mean(ChangeW); CI(ChangeW,prob,2)] ; 
%     DataL=[mean(ChangeL); CI(ChangeL,prob,2)] ; 
%     
%     sumdiff = [sumdiff ; [sum(abs(ChangeW),2) sum(abs(ChangeL),2)]];
%     cat = [cat ; TtB.Treatment(ia)];
%     
%     figure(600)
%     subplot(2,3,t)
%     e=errorbar(Wbins(1:end-1),DataW(1,:),DataW(2,:),DataW(3,:),nanvector,nanvector);
%     e.Color = cmap(t,:);
%     ylabel('Difference between bins [%]')
%     xlabel('Width [mm]')
%     title(Treatments{t})
%     
%     
%     subplot(2,3,3+t)
%     e=errorbar(Lbins(1:end-1),DataL(1,:),DataL(2,:),DataL(3,:),nanvector,nanvector);
%     e.Color = cmap(t,:);
%     ylabel('Difference between bins [%]')
%     xlabel('Length [mm]')
%     title(Treatments{t})
%     
% end
% 
% sgtitle('Difference between Width and length distribution before and after fusion')
% saveas(gcf,'AF_ChangesWLbyBin.png')
% 
% figure;
% subplot(1,2,1)
% boxplot(sumdiff(:,1),cat,'orientation','vertical','label',Treatments,'color',cmap); grid on
% title('Difference of Width distribution')
% xlabel('Treatments'); ylabel('Sum of absolute differences')
% subplot(1,2,2)
% boxplot(sumdiff(:,2),cat,'orientation','vertical','label',Treatments,'color',cmap); grid on
% title('Difference of Length distribution')
% xlabel('Treatments')
% ylabel('Sum of absolute differences')
% saveas(gcf,'AF_TotalWLdifferences.png')

%% 6) Change of Solidity -> Whole cell

% time = 180;
% timesABF = unique(TG.tfusion);
% nanvector = nan(1,numel(timesABF));
% DataSG = zeros(3,numel(timesABF),numel(Treatments)); % Solidity
% 
% diff = [];
% cat ={};
% 
% for t=1:3 % Gather Data
%     Tt = TG(strcmp(TG.Treatment,Treatments{t}),:);    
%     for ati=1:numel(timesABF)
%         Tta = Tt(Tt.tfusion==timesABF(ati),:);
%         S = Tta.At./(Tta.At + Tta.Aws);
%         DataSG(:,ati,t)=[mean(S) CI(S,prob)] ; 
%         
%     end
%     
%     After = Tt(Tt.tfusion==time,:); % After
%     Before = Tt(Tt.tfusion==(-1.*time),:); % Before
%     SA = After.At./(After.At + After.Aws);
%     SB = Before.At./(Before.At + Before.Aws);
%     
%     [~,ia] = intersect([Before.folder Before.Treatment Before.nDish],...
%         [After.folder After.Treatment After.nDish],'rows');
%     
%     diff = [diff ; (SA-SB(ia))./SB(ia)];
%     cat = [cat ; Before.Treatment(ia)];
% 
% end
% 
% figure;
% for t=1:3 % Plot Data
%     subplot(1,2,1)
%     hold on
%     errorbar(timesABF.*(TW/2),DataSG(1,:,t),DataSG(2,:,t),DataSG(3,:,t),nanvector,nanvector)
% end
% 
% subplot(1,2,1)
% ylabel('Solidity [%]')
% xlabel('time from fusion [min]')
% legend(Treatments)
% title('Change in Solidity over time - whole cell')
% 
% subplot(1,2,2)
% boxplot(diff,cat)
% title('Relative Change in Solidity - 4h interval(2 before, 2 after)')
% xlabel('Treatments')
% ylabel('Solidity change [%]')
% 
% saveas(gcf,'AF_GeneralSolidity.png')

%% 7) Change of length and area over time - FusionRegion

% DataA = zeros(3,numel(timesAF),numel(Treatments));
% DataL = zeros(3,numel(timesAF),numel(Treatments));
% 
% A=nan(25,numel(timesAF),3);
% L=nan(25,numel(timesAF),3);
%     
% for t=1:3 % Gather Data
%     Tt = TF(strcmp(TF.Treat,Treatments{t}),:);
%     reps = unique(Tt.Ao);
%       
%     for ati=1:numel(reps)
%         Tta = Tt(Tt.Ao==reps(ati),:);
%         A(ati,:,t) = 100*(Tta.At - Tta.Ao)'/Tta.Ao(1);
%         L(ati,:,t) = 100*(Tta.Lt - Tta.Lo)'/Tta.Ao(1);
%     end
%     DataA(:,:,t)=[nanmean(A(:,:,t),1); CI(A(:,:,t),prob,2)] ; 
%     DataL(:,:,t)=[nanmean(L(:,:,t),1); CI(L(:,:,t),prob,2)] ; 
% end
% 
% idx = [5 9 13];
% figure;
% data = {A(:,idx,1), A(:,idx,2), A(:,idx,3)}; 
% boxplotGroup(data, 'PrimaryLabels', {'C','G','N'}, ...
%   'SecondaryLabels',{'1h after fusion', '2h after fusion','3h after fusion'}, 'GroupLabelType', 'Vertical')
% grid on
% title('Change of Area from fusion [%]')
% saveas(gcf,'AF_AreaChange.png')
% 
% figure;
% data = {L(:,idx,1), L(:,idx,2), L(:,idx,3)}; 
% boxplotGroup(data, 'PrimaryLabels', {'C','G','N'}, ...
%   'SecondaryLabels',{'1h after fusion', '2h after fusion','3h after fusion'}, 'GroupLabelType', 'Vertical')
% grid on
% title('Change of Length from fusion [%]')
% saveas(gcf,'AF_LengthChange.png')


%% % AUXILIAR FUNCTIONS %%%

function [names]=GetSubfolders(path)
        d = dir(path);
        isub = [d(:).isdir]; %# returns logical vector
        names = {d(isub).name}';
        names(ismember(names,{'.','..'})) = [];
end 
function [data] = getFRindexes(Images,treat,folder,times)
    
    nIm = size(Images,3);
    VariableNames = {'Ao' 'Awso' 'Lo' 'At' 'Aws' 'Ac' 'Lt' 'Lv' 'Dt' 'Dv' 'maxD','t'};
    matrix = zeros(nIm,numel(VariableNames));
    
    ppmm =15;
    maxVeinSize = round(1*ppmm)+1; % 1mm meaning 2mm min dimension.
    seE = strel('disk', double(maxVeinSize)); % erode and dilate image
    seD = strel('disk', (round(1.5*maxVeinSize)+1));
    
    for i=1:nIm
        
        I=Images(:,:,i)>0;
        WSA = sum(imfill(I,'holes')-I,'all');
            
        EI = imerode(I,seE);          
        DI = imdilate(EI,seD);

        clustersI = bwlabel(DI.*I); % Get Clusters     
        clusterArea = sum(clustersI,'all');
        
        skel = bwskel(I);
        skel_dist = 2.*(bwdist(~I).*skel);
               
        Ls = [sum(skel,'all') sum(skel.*(~clustersI),'all')];
        Ds = [sum(skel_dist,'all') sum(skel_dist.*(~clustersI),'all')]./Ls;
        
        if i==1 % get initial
            initials = [ sum(I,'all') WSA sum(skel,'all')];
        end
        
        matrix(i,:) = [initials sum(I,'all'), WSA, clusterArea, Ls, Ds, max(skel_dist,[],'all') times(i)];
        
    end
    
    data = array2table(matrix,'VariableNames',VariableNames);
    data = addvars(data,repelem(treat,size(data,1),1),repelem(folder,size(data,1),1),'NewVariableNames',{'Treat','Folder'});
    

end
function [table,ppmm]=getGraphIndexes(Data,times)
  
    nTimes = numel(times); bounds=Data{2};
    IR=zeros(1200,1200); IR(bounds(1):bounds(2),bounds(3):bounds(4)) = Data{3}; IR=IR';
    VariableNames = {'nNodes_tot','nEdges_tot','nNodes_in','nEdges_in','n_Edges_cross','length_path','drag_path','mfn','mfw'};
    var = zeros(nTimes,numel(VariableNames));
    
    for ti=1:nTimes
        
        [G,ppmm] = ConnectGraph(Data{4}{ti});
        
        % Nodes and edges in fusion region (FR)
        nodes_I = accumarray(round([G.Nodes.X,G.Nodes.Y]),1:G.numnodes,[1200 1200]);
        nodes_in = unique(nodes_I.*IR); nodes_in(nodes_in==0)=[];
        edges_in = ismember(G.Edges.EndNodes(:,1),nodes_in) + ismember(G.Edges.EndNodes(:,2),nodes_in);
        
        % Length and drag between centroids
        if ti==1
            c(1,:)=mean([G.Nodes.X(G.Nodes.Type==2) G.Nodes.Y(G.Nodes.Type==2)]);
            c(2,:)=mean([G.Nodes.X(G.Nodes.Type==3) G.Nodes.Y(G.Nodes.Type==3)]);  
            DD = max(max(G.Nodes.Y)-min(G.Nodes.Y),max(G.Nodes.X)-min(G.Nodes.X));
        end
        
        Distance = pdist2(c,[G.Nodes.X G.Nodes.Y]);
        [~,ntc]=min(Distance,[],2);
        Lia = sum(ismember(G.Edges.EndNodes,find((min(Distance,[],1)<DD/4))),2);
        Band = ones(size(Lia)); Band(Lia>1)=100; G.Edges.Weight = Band;
        mfn = maxflow(G,ntc(1),ntc(2));
        Band = G.Edges.Width; Band(Lia>1)=100; G.Edges.Weight = Band;
        mfw = maxflow(G,ntc(1),ntc(2));
        
        G.Edges.Weight = G.Edges.Lp; [~,Dl]=shortestpath(G,ntc(1),ntc(2));
        G.Edges.Weight = G.Edges.drag; [~,Dd]=shortestpath(G,ntc(1),ntc(2));
        
        var(ti,:) = [G.numnodes G.numedges numel(nodes_in) sum(edges_in==2) sum(edges_in==1) Dl Dd mfn mfw];
    end
    
    table = array2table(var,'VariableNames',VariableNames);

end
function [G,m] = ConnectGraph(G)

    % Fix Graph

    Euc_Dist = triu(squareform(pdist([G.Nodes.X G.Nodes.Y])));
    G.Edges.Weight = G.Edges.Le;
    A = triu(full(adjacency(simplify(G),'weighted'))); R = Euc_Dist./A;
    y = A(isfinite(R)); x = Euc_Dist(isfinite(R));
    X = [ones(length(x),1) x]; b = X\y; m=b(2);
    Euc_dist_fixed = b(1) + b(2).*Euc_Dist;
    ifix=find(isnan(G.Edges.Le));
    for j=1:numel(ifix)
        G.Edges.Le(ifix(j)) = Euc_dist_fixed(G.Edges.EndNodes(ifix(j),1),G.Edges.EndNodes(ifix(j),2));         
    end

    G.Edges.Lp(isnan(G.Edges.Lp))=1e-2;
    G.Edges.Lp=max(G.Edges.Lp,G.Edges.Le);

    [bins,binsizes] = conncomp(G);
    if numel(binsizes)>1
    
        o = sortrows([(1:numel(binsizes))' binsizes'],2,"descend");

        edges = zeros(numel(binsizes)-1,3);
        for i=2:numel(binsizes) % Find connecting edge for each subgraph
            gi = o(i,1);
            col_ind = find(bins==gi);
            row_ind = setdiff(1:G.numnodes,col_ind);
            A = pdist2([G.Nodes.X(row_ind) G.Nodes.Y(row_ind)],[G.Nodes.X(col_ind) G.Nodes.Y(col_ind)]);
            [M, I] = min(A(:));
            if M>30
                disp('trouble! more than 30px of distance between subgraphs')
            end
            [r,c]=ind2sub(size(A),I);
            edges(i-1,:) = [row_ind(r) col_ind(c) M*b(2)];  
        end

        NewEdges = table(edges(:,1:2),nanmean(G.Edges.Width).*ones(size(edges,1),1),edges(:,3),edges(:,3),edges(:,3),...
            'VariableNames',{'EndNodes','Width','Le','Lp','Weight'});
        G = addedge(G,NewEdges);

    end
    
    G.Edges.drag = (G.Edges.Lp)./(G.Edges.Width).^4;

    
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
function []=ScatterHistogram(Data,cat,labels,fit)
    % fit = 0 : no fit % fit = 1 : y=mx % fit = 2 : y=mx+c
    
    Treatments = {'Control','Glucose','NaCl'};
    cmap = [27,158,119; 217,95,2; 117,112,179]./255; % Colorblind safe colormap

    figure;
    h = scatterhist(Data(:,1),Data(:,2),'Group',cat,'Marker','+o*','Color',cmap);
    grid on
    xlabel(labels{1})
    ylabel(labels{2})
    
    hold on;
    boxplot(h(2),Data(:,1),cat,'orientation','horizontal','label',Treatments,'color',cmap);
    grid(h(2),'on')
    boxplot(h(3),Data(:,2),cat,'orientation','horizontal','label', Treatments,'color',cmap);
    grid(h(3),'on')
    set(h(2:3),'XTickLabel','');
    view(h(3),[270,90]);  % Rotate the Y plot
    axis(h(1),'auto');  % Sync axes
    
    h(1).Legend.AutoUpdate ="off";
    
    if fit==1
        
        m=Data(:,1)\Data(:,2); rsq=1-sum((Data(:,2)-Data(:,1)*m).^2)/sum((Data(:,2)-mean(Data(:,2))).^2);
        txt1 = ['y = ' num2str(m,3) 'x - Rsq = ' num2str(rsq,2)]; text(0.05, 0.9, txt1,'Units','normalized');
        %text((xL(1)+xL(2))/2,yL(2),txt1,'HorizontalAlignment','left','VerticalAlignment','top','BackgroundColor',[1 1 1],'FontSize',10);
    
    
    elseif fit==2
  
        m=[ones(size(Data,1),1) Data(:,1)]\Data(:,2); rsq=1-sum((Data(:,2)-(Data(:,1)*m(2)+m(1))).^2)/sum((Data(:,2)-mean(Data(:,2))).^2);
        txt1 = ['y = ' num2str(m(2),3) 'x + ' num2str(m(1),3) ' - Rsq = ' num2str(rsq,2)]; text(0.05, 0.9, txt1,'Units','normalized');
          
    end
    
    
    hold off;
    title(labels{3})
end
function[]=plotTreatmentsTime(x,Data,tit,xl,yl,leg)

    cmap = [27,158,119; 217,95,2; 117,112,179]./255;

    figure; hold on

    xr = [x fliplr(x)];
    for t=1:3 % Plot Data
        yr = [Data(1,:,t)+Data(2,:,t), fliplr(Data(1,:,t)+Data(3,:,t))];
        f=fill(xr,yr,cmap(t,:));
        set(f,'facealpha',.25,'LineStyle','none')
        p(t) = plot(x,Data(1,:,t), 'Color',cmap(t,:), 'LineWidth', 2);
    end

    xlabel(xl); ylabel(yl); title(tit); 
    legend([p(1) p(2) p(3)],leg,'Location', 'Best')
    %saveas(gcf,[savepath '.png'])

end

function[]=plotBigraphTime(x,Data1,Data2,tit,xl,yl,leg)
    D{1} = Data1; D{2}= Data2;
    cmap = [27,158,119; 217,95,2; 117,112,179]./255;

    figure;

    xr = [x' fliplr(x')];
    
    for s=1:2
        for t=1:3 % Plot Data
            yr = [D{s}(1,:,t)+D{s}(2,:,t), fliplr(D{s}(1,:,t)+D{s}(3,:,t))];
            subplot(1,2,s); hold on
            f=fill(xr,yr,cmap(t,:));
            set(f,'facealpha',.2,'LineStyle','none')
            p(t) = plot(x,D{s}(1,:,t), 'Color',cmap(t,:), 'LineWidth', 2);
        end
        xlabel(xl); ylabel(yl{s}); 
    end

    sgtitle(tit);
    legend([p(1) p(2) p(3)],leg,'Location', 'Best')
    %saveas(gcf,[savepath '.png'])

end


function[]=plot4Time(x,D,tit,xl,yl,leg)

    cmap = [27,158,119; 217,95,2; 117,112,179]./255;

    figure;

    xr = [x' fliplr(x')];
    
    for s=1:4
        for t=1:3 % Plot Data
            yr = [D{s}(1,:,t)+D{s}(2,:,t), fliplr(D{s}(1,:,t)+D{s}(3,:,t))];
            subplot(2,2,s); hold on
            f=fill(xr,yr,cmap(t,:));
            set(f,'facealpha',.2,'LineStyle','none')
            p(t) = plot(x,D{s}(1,:,t), 'Color',cmap(t,:), 'LineWidth', 2);
        end
        xlabel(xl); ylabel(yl{s}); title(tit{s});
    end

    title(tit);
    legend([p(1) p(2) p(3)],leg,'Location', 'Best')
    %saveas(gcf,[savepath '.png'])

end


function[]=plot3Time(x,D,tit,xl,yl,leg)

    cmap = [27,158,119; 217,95,2; 117,112,179]./255;

    figure;

    xr = [x' fliplr(x')];
    
    for s=1:3
        for t=1:3 % Plot Data
            yr = [D{s}(1,:,t)+D{s}(2,:,t), fliplr(D{s}(1,:,t)+D{s}(3,:,t))];
            subplot(1,3,s); hold on
            f=fill(xr,yr,cmap(t,:));
            set(f,'facealpha',.2,'LineStyle','none')
            p(t) = plot(x,D{s}(1,:,t), 'Color',cmap(t,:), 'LineWidth', 2);
        end
        xlabel(xl); ylabel(yl{s}); title(tit{s});
    end

    title(tit);
    legend([p(1) p(2) p(3)],leg,'Location', 'Best')
    %saveas(gcf,[savepath '.png'])

end