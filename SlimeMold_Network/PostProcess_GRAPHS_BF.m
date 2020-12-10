%clear
close all

% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.065 0.065], [0.2 0.05], [0.075 0.04]); % subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
% if ~make_it_tight,  clear subplot;  end

%load('/Users/lfp3/Dropbox (GaTech)/GT/Spring-20/SMN/Data/A_SampleFolder/G.mat')
%load('/Users/lfp3/Dropbox (GaTech)/GT/Spring-20/SMN/Data/A_SampleFolder/DT.mat')
   
% Constants
cmap = [27,158,119; 217,95,2; 117,112,179]./255; % Colorblind safe colormap
DishTreatments = {'Control','Control','Glucose','Glucose','NaCl','NaCl',...
    'Control','Control','Glucose','Glucose','NaCl','NaCl',...
    'Glucose','Glucose','Glucose','Glucose',...
    'NaCl','NaCl','NaCl','NaCl'};
nDishTreatments = [1 1 2 2 3 3 1 1 2 2 3 3 2 2 2 2 3 3 3 3];

Treatments = {'Control','Glucose','NaCl'};
LegendTreatments = {'Neutral','Nutritive','Adverse'};
shortlegends = {'Neu','Nut','Adv'};

%DataFolder = 'D:\Slime_Mold_Network\Data\Fusion';
%ResultsFolder = 'D:\Slime_Mold_Network\Results\Fusion';

DataFolder = '/Users/lfp3/Dropbox (GaTech)/GT/Spring-20/SMN/Data';
ResultsFolder = '/Users/lfp3/Dropbox (GaTech)/GT/Spring-20/SMN/Data';

% Gather Data
GatherData = false;

if GatherData
    Folders=GetSubfolders(DataFolder);
    Graphs = {{[]} {[]} {[]}};
   
    TableVariables = {'Treatment','Folder','cell','nNodes','nEdges','meanWidth','maxWidth','meanDeg','alpha','tortuosity',...
        'PE_len_m','PE_len_std','PES','PE_drag_m','PE_drag_std',...
        'RR_G_25','RR_G_50','RR_G_75','RR_mst_d_25','RR_mst_d_50','RR_mst_d_75','RR_mst_l_25','RR_mst_l_50','RR_mst_l_75','RR_DT_25','RR_DT_50','RR_DT_75',...
        'RS_wf_25','RS_wf_50','RS_wf_75','RS_sf_25','RS_sf_50','RS_sf_75',...
        'W_all','W_mst_d','W_mst_l','W_DT','BDp_exp','BDp_mean','BLp_exp','BLp_mean'};
    
    varTypes = [{'string'} {'string'} repelem({'double'},numel(TableVariables)-2)];
    sz = [0 numel(TableVariables)];
    T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',TableVariables);
    RData = nan(100,4,150);
    rs = 0;
    for f=1:numel(Folders)
        disp(['Processing folder: :' Folders{f}])
        if isfile(fullfile(DataFolder,Folders{f},'G.mat'))
            S = extractfield(load(fullfile(DataFolder,Folders{f},'G.mat')),'Graphs');
            Gf = S{1};
            for l=1:size(Gf,1)
                t = nDishTreatments(Gf{l,2});
                for c=3:4
                    if ~isempty(Gf{l,c})
                        disp(['l: ' num2str(l) ' cell: ' num2str(c-2)])
                        rs = rs+1;
                        [Data,G,RData(:,:,rs)] = GetGraphData(Gf{l,c},Folders{f},Treatments{t},c-2);
                        if ~isempty(Data)
                            Graphs{t} = [Graphs{t} ;{G}];
                            T = vertcat(T,Data);
                        end
                    end
                end

            end
        end
    end
    
    save(fullfile(ResultsFolder,'GraphsData.mat'),'Graphs','T','RData','-v7.3');
    writetable(T,fullfile(ResultsFolder,'GraphIndexes_BeforeFusion.xlsx'),'Sheet',1,'Range','A1')
else
    load(fullfile(ResultsFolder,'GraphsData.mat'));
end

nReps = cell2mat(cellfun(@numel,Graphs,'UniformOutput',0));

%% START PLOTTING %%

%% Efficiency Decay

is = isnan(permute(RData(1,1,:),[3 1 2]));
RData(:,:,is)=[];
DDD = permute(RData,[3 1 2]);
prctiles = 1:100;
DataRB = zeros(3,numel(prctiles),numel(Treatments));
DataRW = zeros(3,numel(prctiles),numel(Treatments));
DataE = zeros(3,numel(prctiles),numel(Treatments));
DataS = zeros(3,numel(prctiles),numel(Treatments));
prob = 0.95;
for t=1:3 % Gather Data
    inds = strcmp(T.Treatment,Treatments{t});  
    %Tt=T(inds,:); s = Tt.PES;
    for ati=1:numel(prctiles)
        RB = DDD(inds,ati,2);
        DataRB(:,ati,t)=[mean(RB) CI(RB,prob)]; 
        RW = DDD(inds,ati,3);
        DataRW(:,ati,t)=[mean(RW) CI(RW,prob)]; 
        E = DDD(inds,ati,4);%.*s;
        DataE(:,ati,t)=[mean(E) CI(E,prob)]; 
        S = DDD(inds,ati,1);
        DataS(:,ati,t)=[mean(S) CI(S,prob)]; 
    end
    
end

% plotTreatmentsTime(prctiles,DataE,'Efficiency Decay','Percentile [%]','Path Efficiency',Treatments)
% saveas(gcf,'Efficiency_Loss.png'); saveas(gcf,'Efficiency_Loss.fig')


%% GSD 

% figure; hold on
% y=prctiles; Data = DataS;
% yr = [y fliplr(y)];
% for t=1:3 % Plot Data
%     xr = [Data(1,:,t)+Data(2,:,t), fliplr(Data(1,:,t)+Data(3,:,t))];
%     f=fill(xr,yr,cmap(t,:));
%     set(f,'facealpha',.2,'LineStyle','none')
%     p(t) = plot(Data(1,:,t),y, 'Color',cmap(t,:), 'LineWidth', 2);
% end
% 
% xlabel('Drag betweenness'); ylabel('Percentage of edges with lower'); title('Betweeenness Distribution'); 
% legend([p(1) p(2) p(3)],Treatments,'Location', 'Best')
% set(gca, 'XScale', 'log')
% saveas(gcf,'BT_Dist.png'); saveas(gcf,'BT_Dist.fig')

%% Fault Tolerance plot

% Normalized FT
T.NFT_25 = (T.RR_G_25-T.RR_mst_l_25)./(T.RR_DT_25-T.RR_mst_l_25);
T.NFT_50 = (T.RR_G_50-T.RR_mst_l_50)./(T.RR_DT_50-T.RR_mst_l_50);
T.NFT_75 = (T.RR_G_75-T.RR_mst_l_75)./(T.RR_DT_75-T.RR_mst_l_75);
% 
% % Normalized WC
% T.NWCL = (T.W_all-T.W_mst_l)./(T.W_DT-T.W_mst_l);
% 
% 
% cmap = [27,158,119; 217,95,2; 117,112,179]./255;
% figure; hold on
% x=prctiles;
% xr = [x fliplr(x)];
% values = zeros(3,3);
% for t=1:3 % Plot Data
%     yr = [DataRB(1,:,t)+DataRB(2,:,t), fliplr(DataRB(1,:,t)+DataRB(3,:,t))];
%     f=fill(xr,yr,cmap(t,:));
%     set(f,'facealpha',.2,'LineStyle','none')
%     p(t) = plot(x,DataRB(1,:,t), 'Color',cmap(t,:), 'LineWidth', 2);
%     
%     yr = [DataRW(1,:,t)+DataRW(2,:,t), fliplr(DataRW(1,:,t)+DataRW(3,:,t))];
%     f=fill(xr,yr,cmap(t,:));
%     set(f,'facealpha',.2,'LineStyle','none')
%     plot(x,DataRW(1,:,t), 'Color',cmap(t,:), 'LineWidth', 2);
%     
%     values(:,t) = [100*mean(T.NFT_25(strcmp(T.Treatment,Treatments{t})));...
%         100*mean(T.NFT_50(strcmp(T.Treatment,Treatments{t})));...
%         100*mean(T.NFT_75(strcmp(T.Treatment,Treatments{t})))];
%     
%     scatter(values(1,t),25,[],cmap(t,:),'filled');
%     scatter(values(2,t),50,[],cmap(t,:),'filled');
%     scatter(values(3,t),75,[],cmap(t,:),'filled');
%     
% end
% 
% yval = 0:25:100;
% values = [100 100 100; values; 0 0 0];
% for t=1:3
%     yy = interp1(values(:,t),yval,0:100,'pchip');
%     l = plot(0:100,yy,'-.','Color',cmap(t,:));
% end
% 
% plot([0 100],[25 25],':k'); plot([0 100],[50 50],':k'); plot([0 100],[75 75],':k')
% xlim([0 100]); ylim([0 100])
% xlabel('Removed Edges [%]'); ylabel('Connected Nodes [%]'); title('Fault tolerance by Edge Importance'); 
% legend([p(1) p(2) p(3)],LegendTreatments,'Location', 'Best')
% 
% axxx=axes('position',get(gca,'position'),'visible','off');
% legend(axxx,l,'Mean Random FT','Location','EastOutside');
% 
% text(0.55, 0.6,'Important Last','Units','normalized');
% text(0.05, 0.3,'Important First','Units','normalized');
% saveas(gcf,'FT_Plot.png'); saveas(gcf,'FT_Plot.fig')
% print('FT_Plot','-depsc');

%% Robustness / Fault Tolerance

% figure; 
% subplot(1,3,1);
% boxplot(T.NFT_75,T.Treatment,'orientation','vertical','label',shortlegends,'color',cmap); 
% ylabel('NFT - 75% Connected'); grid on; ylim([0 0.5])
% subplot(1,3,2);
% boxplot(T.NFT_50,T.Treatment,'orientation','vertical','label',shortlegends,'color',cmap); 
% ylabel('NFT - 50% Connected'); grid on; ylim([0 0.5])
% subplot(1,3,3);
% boxplot(T.NFT_25,T.Treatment,'orientation','vertical','label',shortlegends,'color',cmap); 
% ylabel('NFT - 25% Connected'); grid on; ylim([0 0.5])
% 
% sgtitle('Normalized Fault tolerance - NFT [-]')
% 
% saveas(gcf,'NormalizedFaultTolerances.png'); saveas(gcf,'NormalizedFaultTolerances.fig')
% print('NormalizedFaultTolerances','-depsc');

% RD = [T.RS_sf./T.RR_G, T.RS_wf./T.RR_G]; 
% data = {RD(strcmp(T.Treatment,Treatments{1}),:),RD(strcmp(T.Treatment,Treatments{2}),:),RD(strcmp(T.Treatment,Treatments{3}),:)};
% subplot(1,2,2); boxplotGroup(data, 'PrimaryLabels', {'C' 'G' 'S'}, ...
%   'SecondaryLabels',{'Important Edges First','Important Edges Last'},...
%   'GroupLabelType', 'Vertical'); ylabel('Hierarchical Fault Tolerance [-]'); grid on
% saveas(gcf,'BP_FT.png'); saveas(gcf,'BP_FT.fig')


% % figure; subplot(1,2,1)
% % y=T.NFT_50; x=T.NWCL;
% % b=[ones(numel(x),1) x]\y; rsq=1-sum((y-(x*b(2)+b(1))).^2)/sum((y-mean(y)).^2);
% % hold on
% % txt1 = sprintf(['NFT = ' num2str(b(2),3) ' NNL + ' num2str(b(1),3) ' \n R^2 = ' num2str(rsq,2)]); text(0.05, 0.1, txt1,'Units','normalized');
% % plot([0 max(x)],b(1)+b(2).*[0 max(x)],':k')
% % gscatter(x,y,T.Treatment,cmap,'*+o')
% % xlabel('Normalized Network Length (NNL)'); ylabel('Normalized Fault tolerance (NFT)')
% % title('Network Cost vs Resiliency')
% % 
% % subplot(1,2,2)
% % y=T.BDp_exp; x=T.NWCL; %y=T.PES ; %
% % b=[ones(numel(x),1) x]\y; rsq=1-sum((y-(x*b(2)+b(1))).^2)/sum((y-mean(y)).^2);
% % hold on
% % txt1 = sprintf(['MEB = ' num2str(b(2),3) ' NNL + ' num2str(b(1),3) ' \n R^2 = ' num2str(rsq,2)]); text(0.05, 0.1, txt1,'Units','normalized');
% % plot([0 max(x)],b(1)+b(2).*[0 max(x)],':k')
% % gscatter(x,y,T.Treatment,cmap,'*+o'); ylim([0 0.1])
% % ylabel('Mean Edge Betweenness (MEB)'); xlabel('Normalized Network Length (NNL)')
% % title('Network Cost vs. Edge Betweenness');
%axis equal


% subplot(1,2,2);
% boxplot(y./x,T.Treatment,'orientation','vertical','label',Treatments,'color',cmap); grid on
% title('Path efficiency over Fault Tolerance')

% saveas(gcf,'Relationship_PE_FT.png'); saveas(gcf,'Relationship_PE_FT.fig')
% print('Relationship_PE_FT','-depsc');

% Robustness range
% RR = T.RS_sf./T.RS_wf;
% figure;
% boxplot(RR,T.Treatment)


% figure; hold on
% for t=1:3
%     inds = strcmp(T.Treatment,Treatments{t});  
%     A = RData(:,:,inds);
%     C = permute(A,[1 3 2]); C = reshape(C,[],size(A,2),1); C = sortrows(C,1);    
%     fitresult = fit(C(:,1),C(:,2),'exp1'); p=plot(fitresult); set(p,'color',cmap(t,:))
%     %p22 = predint(fitresult,C(:,1),0.95,'functional','off'); plot(C(:,1),p22,'m--')
%     fitresult = fit(C(:,1),C(:,3),'exp1');  p=plot(fitresult); set(p,'color',cmap(t,:))
%     %p22 = predint(fitresult,C(:,1),0.95,'functional','off'); plot(C(:,1),p22,'m--')    
% end

% figure; hold on
% for t=1:3
%     inds = strcmp(T.Treatment,Treatments{t});  
%     A = RData(:,:,inds);
%     C = permute(A,[1 3 2]); C = reshape(C,[],size(A,2),1); C = sortrows(C,1);    
%     [fitresult,gof,~] = fit(C(:,1),C(:,4),'exp1'); p=plot(fitresult); set(p,'color',cmap(t,:))
%     gof.rsquare  
% end

%% nodes vs edges
% % ScatterHistogram([T.nNodes T.nEdges],T.Treatment,{'Number of Nodes','Number of Edges',''},1);
% % hold on; plot([0 1000],[0 1000].*(T.nNodes\T.nEdges),':k')
%saveas(gcf,'SH_numNodes_numEdges.png'); saveas(gcf,'SH_numNodes_numEdges.fig')
%print('SH_numNodes_numEdges','-depsc');

% Alpha - Degree
figure; subplot(1,2,1); boxplot(T.alpha,T.Treatment,'orientation','vertical','label',Treatments,'color',cmap); ylabel('Alpha'); grid on
subplot(1,2,2); boxplot(T.meanDeg,T.Treatment,'orientation','vertical','label',Treatments,'color',cmap); ylabel('Mean node Degree'); grid on
sgtitle('Graph connectivity: Alpha and Node degree')
% saveas(gcf,'BP_Connectivity.png'); saveas(gcf,'BP_Connectivity.fig')


%% Wiring Cost - Tortuosity
% ScatterHistogram([T.W_mst_l/10 T.W_mst_d/10],T.Treatment,{'MST - Minimum Length [cm]','MST - Minimum Drag [cm]','MST Network length'},1);
% hold on; plot([0 100],[0 100].*(T.W_mst_l\T.W_mst_d),':k')
% %legend([Treatments {['m=' num2str(round(T.W_mst_l\T.W_mst_d,2),3)]}])
% saveas(gcf,'SH_MST_lengths.png'); saveas(gcf,'SH_MST_lengths.fig')


T.Normalized_Network_Length = (T.W_all-T.W_mst_l)./(T.W_DT-T.W_mst_l);
% Ratio - Tortuosity
figure; subplot(1,2,1); boxplot(T.tortuosity,T.Treatment,'orientation','vertical','label',LegendTreatments,'color',cmap); 
ylabel('Tortuosity [%]'); ylim([0 6]); grid on
title('a) Network tortuosity')
subplot(1,2,2); boxplot(T.Normalized_Network_Length,T.Treatment,'orientation','vertical','label',LegendTreatments,'color',cmap); 
ylabel('Normalized network length [-]'); grid on; ylim([0.075 0.225])
title('b) Normalized network length')
% saveas(gcf,'BP_Tortuosity-Length.png'); saveas(gcf,'BP_Tortuosity-Length.fig')
% print('BP_Tortuosity-Length','-depsc');

%% Path Efficiency
T.PE_len_cov = T.PE_len_std./T.PE_len_m;
[filtered,tf] = rmoutliers([T.PE_len_std./T.PE_len_m T.PE_len_m]);
% ScatterHistogram(filtered,T.Treatment(~tf),{'Coefficient of Variation (COV) [-]','Mean [-]','LE: Length efficiency distribution'},0);
% saveas(gcf,'SH_LE.png'); saveas(gcf,'SH_LE.fig')

figure; subplot(1,4,1); boxplot(filtered(:,2),T.Treatment(~tf),'orientation','vertical','label',shortlegends,'color',cmap); 
ylabel('Mean [-]'); grid on
subplot(1,4,2); boxplot(filtered(:,1),T.Treatment(~tf),'orientation','vertical','label',shortlegends,'color',cmap); 
ylabel('Coefficient of variation (COV) [-]'); grid on
title('a) Length efficiency distribution (LE)')
% saveas(gcf,'BP_LE.png'); saveas(gcf,'BP_LE.fig')
% print('BP_LE','-depsc');

T.PE_drag_cov = T.PE_drag_std./T.PE_drag_m;
[filtered,tf] = rmoutliers([T.PE_drag_std./T.PE_drag_m T.PE_drag_m]);
% ScatterHistogram(filtered,T.Treatment(~tf),{'Coefficient of Variation (COV) [-]','Mean [-]','DE: Drag Efficiency distribution'},0);
% saveas(gcf,'SH_DE.png'); saveas(gcf,'SH_DE.fig')

subplot(1,4,3); boxplot(filtered(:,2),T.Treatment(~tf),'orientation','vertical','label',shortlegends,'color',cmap); 
ylabel('Mean [-]'); grid on
subplot(1,4,4); boxplot(filtered(:,1),T.Treatment(~tf),'orientation','vertical','label',shortlegends,'color',cmap); 
ylabel('Coefficient of variation (COV) [-]'); grid on
title('b) Drag efficiency distribution (DE)')
% saveas(gcf,'BP_DE.png'); saveas(gcf,'BP_DE.fig')
% print('BP_DE','-depsc');

% LE vs DE
% ScatterHistogram([T.PE_len_m T.PE_drag_m],T.Treatment,{'LE: Length efficiency','DE: Drag Efficiency','Drag vs Length Efficiency'},0);
% xlim([0.5 0.9]);  ylim([0 15])
% saveas(gcf,'SH_LEvsDE.png'); saveas(gcf,'SH_LEvsDE.fig')

%% Betweenness

figure; subplot(1,2,1); boxplot(T.BLp_exp,T.Treatment,'orientation','vertical','label',LegendTreatments,'color',cmap); 
ylabel('BL: Mean edge betweenness [%]'); grid on; ylim([0 0.12]); title('Betweenness by length (BL)')
subplot(1,2,2); boxplot(T.BDp_exp,T.Treatment,'orientation','vertical','label',LegendTreatments,'color',cmap); 
ylabel('BD: Mean edge betweenness [%]'); grid on; ylim([0 0.12]); title('Betweenness by drag (BD)')
% 
% saveas(gcf,'BP_BD_BE.png'); saveas(gcf,'BP_BD_BE.fig')
% print('BP_BD_BE','-depsc');

% [filtered,tf] = rmoutliers([T.BDp_exp T.BLp_exp]);
% ScatterHistogram(filtered,T.Treatment(~tf),{'BD: Mean edge Betweenness by drag [%]','BL: Mean edge Betweenness by length [%]','Edge Betweenness distribution'},2);
% hold on; plot([0 0.1],[0 0.1].*(filtered(:,1)\filtered(:,2)),':k')
% saveas(gcf,'SH_BD_BE.png'); saveas(gcf,'SH_BD_BE.fig')

% figure;
% x=(T.W_all-T.W_mst_l)./(T.W_DT-T.W_mst_l); y=T.PE_len_m;
% b=[ones(numel(x),1) x]\y; rsq=1-sum((y-(x*b(2)+b(1))).^2)/sum((y-mean(y)).^2);
% hold on
% txt1 = ['y = ' num2str(b(2),3) 'x + ' num2str(b(1),3) ' - Rsq = ' num2str(rsq,2)]; text(0.05, 0.9, txt1,'Units','normalized');
% plot([0 max(x)],b(1)+b(2).*[0 max(x)],':k')
% gscatter(x,y,T.Treatment,cmap,'*+o')
% ylabel('Path Length Efficiency'); xlabel('Normalized Network Length')
% title('Network Length vs Path Efficiency')

%figure; boxplot(y./x,T.Treatment); title('Path efficiency over network length')

save(fullfile(ResultsFolder,'GraphsData.mat'),'Graphs','T','RData','-v7.3');
writetable(T,fullfile(ResultsFolder,'GraphIndexes_BeforeFusion.xlsx'),'Sheet',1,'Range','A1')

%% AUXILIAR FUNCTIONS %%

function [names]=GetSubfolders(path)
        d = dir(path);
        isub = [d(:).isdir]; %# returns logical vector
        names = {d(isub).name}';
        names(ismember(names,{'.','..'})) = [];
end 
function []=plotGraph(G,edgePropToColor)
    %map = brewermap(101,'YlOrRd');
    map = copper(101);
    scaledV= round(100*edgePropToColor./prctile(edgePropToColor,95));
    scaledV(scaledV>100)=100;
    
    NodeC=brewermap(3,'Set1');
    Colors = map(scaledV+1,:);
    figure;
    %plot(G,'EdgeColor',Colors,'NodeColor',NodeC(G.Nodes.Type+1,:),'NodeLabel',{});
    plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'EdgeColor',Colors,'NodeColor',NodeC(G.Nodes.Type+1,:),'NodeLabel',{});
    
    xlabel('X Coordinate')
    ylabel('Y Coordinate')
    axis equal
    axis tight

end
function [W] = WiringCost(G)
    G.Edges.Weight = G.Edges.Lp;
    adj = adjacency(simplify(G),'weighted');
    W = sum(adj(:))/2;
end
function [T,G,Data] = GetGraphData(G,Folder,Treatment,cell)

    G = simplify(G);
    if G.numedges<10
        T=[];
        Data=nan(100,4);
        disp('Less than 10 edges! NOT COMPUTED');
        return
    end
    % Join disconnected graph
    [~,binsizes] = conncomp(G);   
    if numel(binsizes)>1; G = ConnectGraph(G); end
    
    % Fix Coordinates and lengths
    
    G.Edges.Le(isnan(G.Edges.Le))=1e-1;
    G.Edges.Lp(isnan(G.Edges.Lp))=1e-1;
    
    % Fix Coordinates    
    Euc_Dist = triu(squareform(pdist([G.Nodes.X G.Nodes.Y])));
    G.Edges.Weight = G.Edges.Le;
    A = triu(full(adjacency(simplify(G),'weighted')));
    R = Euc_Dist./A;
    y = A(isfinite(R));
    x = Euc_Dist(isfinite(R));
    X = [ones(length(x),1) x];
    b = X\y;
    %figure; scatter(Euc_Dist(isfinite(R)),A(isfinite(R)))
    G.Nodes.X = b(2).*G.Nodes.X; G.Nodes.Y = b(2).*G.Nodes.Y;
    G.Edges.Le = G.Edges.Le - b(1);
    G.Edges.Lp=max(G.Edges.Lp,G.Edges.Le);
    G.Edges.drag = (G.Edges.Lp)./(G.Edges.Width).^4;
    
    %% Betweenness and path efficiency
    %AP=find(or(G.Nodes.Type==0,G.Nodes.Type>2));
    %AP=find(G.Nodes.Type>0);
    AP = 1:G.numnodes;

    GD = G; 
    G.Edges.Weight = G.Edges.Lp; PP = triu(distances(G));
    GD.Edges.Weight = G.Edges.drag; PD = triu(distances(GD));
    
    edge_count_D = zeros(size(G.Edges,1),1);
    edge_count_P = zeros(size(G.Edges,1),1);
    
    % Loop computing shortest paths
    for s=1:length(AP)
        for t=s+1:length(AP)
           [~,~,eu] = shortestpath(GD,AP(s),AP(t)); 
           edge_count_D(eu) = edge_count_D(eu)+1;
           [~,~,eu] = shortestpath(G,AP(s),AP(t)); 
           edge_count_P(eu) = edge_count_P(eu)+1;
        end
    end

    G.Edges.BL = (2./((numel(AP)-2)*(numel(AP)-1))).*edge_count_D;
    G.Edges.BD = (2./((numel(AP)-2)*(numel(AP)-1))).*edge_count_P;
     
    Bd = fitdist(G.Edges.BD,'exponential'); BDp = [Bd.mu Bd.mean];
    Bd = fitdist(G.Edges.BL,'exponential'); BLp = [Bd.mu Bd.mean];
    
    %plotGraph(G,G.Edges.BL)
    %plotGraph(G,G.Edges.BD)
    
    % BN
%     n = numnodes(G);
%     wbc = centrality(G,'betweenness','Cost',G.Edges.Weight+0.0001);
%     NC = 2*wbc./((n-2)*(n-1));
%     G.Edges.BN = min(NC(G.Edges.EndNodes),[],2);


    %% Compute path efficiency: Drag and Length
    Euc_Dist = triu(squareform(pdist([G.Nodes.X G.Nodes.Y])));
    ELM = Euc_Dist./PP; % Efficiency Length
    EL = ELM(isfinite(ELM)); EL(EL>1)=1;
    %figure; histogram(EL)
    E=PathEfficiencyScalar(G);
    EL = [mean(EL) std(EL) E];    
    DV = 0.5^(-4).*(Euc_Dist); DV(DV<=0)=[]; %prctile(G.Edges.Width,50)
    PDV = PD(PD>0);
    EDM = DV./PDV'; EDM(EDM>prctile(EDM,95))=[];
    ED = [mean(EDM) std(EDM)];
    %edges = linspace(0,prctile(EDM,95),30); figure; histogram(EDM,edges)
 
    %% Wiring Cost
    
    W = [0 0 0];
    W(1) = WiringCost(G);
    
    G.Edges.Weight = G.Edges.drag;
    MST_d = minspantree(G);
    W(2) = WiringCost(MST_d);
    
    G.Edges.Weight = G.Edges.Lp;
    MST_l = minspantree(G);
    W(3) = WiringCost(MST_l);
       
    tri = delaunay([G.Nodes.X G.Nodes.Y]);
    A = adjacency(digraph(tri, tri(:, [2 3 1]))); A = A | A';
    w = A.*(Euc_Dist+Euc_Dist');
    DT = graph(w); DT.Edges.Lp = DT.Edges.Weight;
    W(4) = WiringCost(DT);
        
%     AP = boundary(G.Nodes.X,G.Nodes.Y,0.5);
%     %AP = find(G.Nodes.Type>0); 
%     counts = zeros(1,MST.numnodes);
%     for s=1:length(AP)
%         for t=s+1:length(AP)
%            P = shortestpath(MST,AP(s),AP(t)); 
%            if isempty(P); continue; end
%            counts = counts + histc(P,1:MST.numnodes);
%         end
%     end
%     MST = rmnode(MST,find(counts==0));

    %% ROBUSTNESS
    
    % Random
    RR(1:3)=RobustRandom(G);
    RR(4:6)=RobustRandom(MST_d);
    RR(7:9)=RobustRandom(MST_l);
    RR(10:12)=RobustRandom(DT);
    
    % Sorted - Remove Weak
    metric = G.Edges.BD;
    Data = [prctile(metric,1:1:100)' zeros(100,3)] ;
    for i=1:size(Data,1)
        idx = find(metric<=Data(i,1));
        [~,binsize] = conncomp(rmedge(G,idx));
        Data(i,2) = 100*(max(binsize)/G.numnodes);
    end 
    
    RS(1:3) = [find(Data(:,2)<25,1) find(Data(:,2)<50,1) find(Data(:,2)<75,1)] ;
    
    %figure; plot(1:100,Data(:,2))
    
    % Sorted - Remove Strong
    for i=1:size(Data,1)
        idx = find(metric>Data(size(Data,1)-i+1,1));
        GR = rmedge(G,idx);
        Data(i,4) = PathEfficiencyScalar(GR);
        [~,binsize] = conncomp(GR);
        Data(i,3) = 100*(max(binsize)/G.numnodes);
    end 
    Data(:,4)=Data(:,4)./E;
    RS(4:6) = [find(Data(:,3)<25,1) find(Data(:,3)<50,1) find(Data(:,3)<75,1)] ;
    
    %hold on; plot(1:100,Data2(:,2))
    %plot([0 100],[50 50])
    
    % Mean-Max Vein Width 
    meanV = nanmean(G.Edges.Width);
    maxV = nanmax(G.Edges.Width);
    
    %% Fill data and return
    ND = [cell,G.numnodes,G.numedges,meanV,maxV,mean(degree(G)),(G.numedges-G.numnodes+1)/(2*G.numedges-5),...
        sum(G.Edges.Lp)/sum(G.Edges.Le),EL,ED,RR,RS,W,BDp,BLp];
    % [Cell, #Nodes, #Edges, meanDegree, alpha,tortuosity,
    % EfLen_m,Eflen_std,EfDrag_m,EfDrag_std,
    % RobustRandom,RobustSort_wf,RobustSort_sf,W_all,W_mst,W_tt,BDp_exp,BDp_mean,BLp_exp,BLp_mean]
    
    TableVariables = {'cell','nNodes','nEdges','meanWidth','maxWidth','meanDeg','alpha','tortuosity',...
        'PE_len_m','PE_len_std','PES','PE_drag_m','PE_drag_std',...
        'RR_G_25','RR_G_50','RR_G_75','RR_mst_d_25','RR_mst_d_50','RR_mst_d_75','RR_mst_l_25','RR_mst_l_50','RR_mst_l_75','RR_DT_25','RR_DT_50','RR_DT_75',...
        'RS_wf_25','RS_wf_50','RS_wf_75','RS_sf_25','RS_sf_50','RS_sf_75',...
        'W_all','W_mst_d','W_mst_l','W_DT','BDp_exp','BDp_mean','BLp_exp','BLp_mean'};

    T = array2table(ND,'VariableNames',TableVariables);
    T = addvars(T,{Treatment},{Folder},'NewVariableNames',{'Treatment','Folder'});

end
function [G] = ConnectGraph(G)

    [bins,binsizes] = conncomp(G);
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
        edges(i-1,:) = [row_ind(r) col_ind(c) M];  
    end
    
    NewEdges = table(edges(:,1:2),nanmean(G.Edges.Width).*ones(size(edges,1),1),edges(:,3),edges(:,3),...
        'VariableNames',{'EndNodes','Width','Le','Lp'});
    G = addedge(G,NewEdges);

end

function [E]=PathEfficiencyScalar(G)
    PP = 1./triu(distances(G));
    Euc_Dist = 1./triu(squareform(pdist([G.Nodes.X G.Nodes.Y])));
    E = sum(PP(isfinite(PP)))/sum(Euc_Dist(isfinite(PP))); % Efficiency Length    
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
    set(f,'facealpha',.2,'LineStyle','none')
    p(t) = plot(x,Data(1,:,t), 'Color',cmap(t,:), 'LineWidth', 2);
end

xlabel(xl); ylabel(yl); title(tit); 
legend([p(1) p(2) p(3)],leg,'Location', 'Best')
%saveas(gcf,[savepath '.png'])

end


function [RR]=RobustRandom(G)

    nReps = 30;
    R=zeros(nReps,1);
    
    for r=1:nReps
        idxs = randperm(G.numedges);
        i=0;
        pE=100;
        while pE>25
            i=i+1;
            idx = idxs(1:round(G.numedges*i/100));
            [~,binsize] = conncomp(rmedge(G,idx));
            pE = 100*(max(binsize)/G.numnodes);
        end 
        R(r) = i;
    end
    
    RR(1)=nanmean(R);
    

    R=zeros(nReps,1);
    
    for r=1:nReps
        idxs = randperm(G.numedges);
        i=0;
        pE=100;
        while pE>50
            i=i+1;
            idx = idxs(1:round(G.numedges*i/100));
            [~,binsize] = conncomp(rmedge(G,idx));
            pE = 100*(max(binsize)/G.numnodes);
        end 
        R(r) = i;
    end
    
    RR(2)=nanmean(R);
    
     R=zeros(nReps,1);
    
    for r=1:nReps
        idxs = randperm(G.numedges);
        i=0;
        pE=100;
        while pE>75
            i=i+1;
            idx = idxs(1:round(G.numedges*i/100));
            [~,binsize] = conncomp(rmedge(G,idx));
            pE = 100*(max(binsize)/G.numnodes);
        end 
        R(r) = i;
    end
    
    RR(3)=nanmean(R);
      
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

