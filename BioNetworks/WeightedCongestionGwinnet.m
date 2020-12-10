% Weighted congestion

path='/Users/lfp3/Dropbox (GaTech)/GT/Spring-19/RootsPaper/ATL_map_files/Gwinett_WT_Dens_Centroids_N5.csv';
%path  = 'C:\Users\lfp3\Dropbox (GaTech)\GT\Spring-19\RootsPaper\ATL_map_files\Gwinett_WT_Dens_Centroids_N5.csv';
%Points=dlmread(path,',',1,0);
Points(1,3)=3639843-sum(Points(2:end,3));
%Points(1,3)=sum(Points(2:end,3));

a=Points(:,1:2); % a has to have just coordinates
points=a;

% Fan Network
gg=load('Gwinnett_FanGraph_N5.mat');
GF=gg.GF;
% Steiner Network:
gg=load('Gwinnett_SteinerGraph_N5.mat');
GS=gg.GS;
% Leaf Network
gg=load('Gwinnett_LeafGraph_N5_2.mat');
GL=gg.GL;

%ojjhrder = [7 8 2 9 5 1 3 4 10 6];
%order = [6 9 7 8 5 10 1 2 4 3];
%GL = reordernodes(GL,order);

Gs={GF,GS,GL};
Labels={'Fan','Steiner','Leaf'};

weights = zeros(6,1);
OrderedPoints = zeros(6,2);
EC = cell(2,1);
Dist = cell(2,1);

% Row 1 - ST / Row 2 - LV
%Cols: Population, Area, Length

GIS=[231496 199.38890609 47.185;254594 207.81981171 50.052];


for i=2:3
    G=Gs{i};
    
    AP=find(or(G.Nodes.Type==1,G.Nodes.Type==3));
    pG = [G.Nodes.X G.Nodes.Y];
    
    for k=1:numel(AP)      
        [v,I]=min(pdist2(points,pG(AP(k),:)));
        I
        weights(k)=Points(I,3);
        OrderedPoints(k,:)=Points(I,1:2);
    end
    edges_used=[];

    for s=1:length(AP)
        for t=s+1:length(AP)
           [P,d,eu] = shortestpath(G,AP(s),AP(t)); 
           ew = [eu' repmat((weights(s)+weights(t))/2,numel(eu),1)];
           edges_used=[edges_used; ew];
        end
    end
    
    edge_count = zeros(size(G.Edges,1),1);
        
    for e=1:size(G.Edges,1)
        edge_count(e) = sum(edges_used(edges_used(:,1)==e,2));
    end

    EC{i-1} = edge_count./sum(weights);

    
    Euc_Dist = squareform(pdist(OrderedPoints)); % Euclidean Distance between pairs of points

    % PairWise Distances
    d=triu(distances(G,AP,AP)./Euc_Dist); %km
    d(d==0)=nan;
    
    Dist{i-1}=d;

end

%% Dist Ratio

figure(134)
subplot(1,2,1)
boxplot([Dist{2}(1,:)' Dist{1}(1,:)'],'Labels',{'LV','ST'})
title('Path ratio - From source')
grid on
subplot(1,2,2)
boxplot([Dist{2}(:) Dist{1}(:)],'Labels',{'LV','ST'})
title('Path ratio - All Centroids')
grid on

disp('All to all')
disp('Steiner: mean - iqr')
[nanmean(Dist{1}(:)) iqr(Dist{1}(:))]
disp('Leaf: mean - iqr')
[nanmean(Dist{2}(:)) iqr(Dist{2}(:))]

disp('From Center')
disp('Steiner: mean - iqr')
[nanmean(Dist{1}(1,:)) iqr(Dist{1}(1,:))]
disp('Leaf: mean - iqr')
[nanmean(Dist{2}(1,:)) iqr(Dist{2}(1,:))]


%% Weighted Congestion:

mean(EC{1})
mean(EC{2})
iqr(EC{1})
iqr(EC{2})

figure(258)
subplot(1,2,2)

pe = plot(GS,'XData',GS.Nodes.X,'YData',GS.Nodes.Y,'NodeLabel',{},'EdgeAlpha',1,'LineWidth',EC{1});
pe.EdgeCData = EC{1};
pe.MarkerSize=1;
colormap(flip(winter,1));
set(gca,'XTick',[]);
set(gca,'YTick',[]);
title('Steiner Tree')

subplot(1,2,1)

pe = plot(GL,'XData',GL.Nodes.X,'YData',GL.Nodes.Y,'NodeLabel',{},'EdgeAlpha',1,'LineWidth',EC{2});
pe.EdgeCData = EC{2};
pe.MarkerSize=1;
colormap(flip(winter,1));
set(gca,'XTick',[]);
set(gca,'YTick',[]);
title('Leaf Venation')




function [I]=GetIndexes(G)
%{
1. Total Length - Sum dis_from_center
2. Mean ratio Path/Euclidean_dist -> From Center
3. Mean ratio Path/Euclidean_dist -> From every point
4.  Node Degree
5.  Betweeness
6.  Congestion
7.  Page rank
%}

    I=zeros(2,2);
    pi = find(~(G.Nodes.Type==2));
    points = [G.Nodes.X(pi) G.Nodes.Y(pi)];

    Euc_Dist = squareform(pdist(points)); % Euclidean Distance between pairs of points

    % PairWise Distances
    d=triu(distances(G,pi,pi)./Euc_Dist); %km
    d(d==0)=nan;
    
    I(1,1)=sum(G.Edges{:,2}); %Total network length
    I(1,2)=sum(distances(G,1,pi)); %Sum of distances from the source
    I(2,:)=[nanmean(d(1,:)) iqr(d(1,:))]; % Sum of
    I(3,:)=[nanmean(d(:)) iqr(d(:))];
    I(4,:)=[mean(degree(G)) iqr(degree(G))];
    I(5:end,:)=GraphImportance(G);
    
end