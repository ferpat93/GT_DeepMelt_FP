% Graph to Lines QGIS

% Format: [X Y Weight Order Group]

% Fan Network
gg=load('FanGraph_N15.mat');
GF=gg.GF;
% Steiner Network:
gg=load('SteinerGraph_N15.mat');
GS=gg.GS;
% Leaf Network
gg=load('LeafGraph_N15.mat');
GL=gg.GL;

Gs={GF,GS,GL};
Labels={'Fan','Steiner','Leaf'};

%{
write_lines=0;

if write_lines==1   
    for i=1:numel(Gs)

        path=fullfile(pwd,'ATL_map_files',strcat('Graph_Line_Layer_',Labels{i},'.csv'));

        G=Gs{i};
        Nodes=G.Nodes{:,:};
        Edges=G.Edges{:,:};

        X = reshape([Nodes(Edges(:,1),1)';Nodes(Edges(:,2),1)'], 1, [])';
        Y = reshape([Nodes(Edges(:,1),2)';Nodes(Edges(:,2),2)'], 1, [])';
        order=repmat([1;2],size(Edges,1),1);
        group = reshape([1:size(Edges,1);1:size(Edges,1)], 1, [])';

        T = table(X,Y,order,group);
        writetable(T,path,'WriteRowNames',true) 

    end
end
%}

%% Get Data
% Order: 1.Fan  2.Steiner   3.Leaf  4.Roots   5.Marta -> Columns
    % Fan: ﻿1228947 ﻿940552893.1543258
    % Steiner: ﻿857522 ﻿650758270.4742482
    % leaf: ﻿935496 ﻿712078782.340021
    % Marta: ﻿479723 ﻿295356142
    
% Rows: 
%{
1. Population
2. Buffer Area (2km)
3. Total Length
4. Sum dis_from_center
5. Buffer Density (pop/area)
6. Efficiency(?): (pop/total_length)
7. Mean ratio Path/Euclidean_dist -> From Center
8. Mean ratio Path/Euclidean_dist -> From every point
9. Buff Area/Total Length
%}

Data=zeros(10,4);
Data(1:2,5)=[479723 295356142]';
Data(1:2,1)=[1228947 940552893.1543258]';
Data(1:2,3)=[935496 712078782.340021]';
Data(1:2,2)=[857522 650758270.4742482]';

Data(2,:)=Data(2,:)./1000000;

% Coordinates of 15 attraction points
points=GF.Nodes{:,1:2};

Euc_Dist = squareform(pdist(points)); % Euclidean Distance between pairs of points
%triu()

% PairWise Distances
pw_dist=zeros(16,16,3);

%Importance data
IData=zeros(3,4);

for i=1:3
    G=Gs{i};
    d=triu(distances(G,1:16,1:16)./Euc_Dist)./1000; %km
    d(d==0)=nan;
    pw_dist(:,:,i)=d;
    
    Data(3,i)=sum(G.Edges{:,2})/1000; %km
    Data(4,i)=sum(distances(G,1,2:16))/1000;%km
    Data(5,i)=Data(1,i)/Data(2,i);
    Data(6,i)=Data(1,i)/Data(3,i);
    Data(7,i)=nanmean(d(1,:));
    Data(8,i)=nanmean(d(:));
    Data(9,i)=Data(2,i)/Data(3,i);
    
    %{
    figure(i+10)
    subplot(1,2,1);  
    heatmap(d,'Colormap',summer);
    title('Pairwise HeatMap')
    subplot(1,2,2);
    histogram(d,'FaceColor','r');
    title('Frequency Histogram')   
    %sgtitle([Labels{i} ' network: Ratio between travel distance and euclidean distance'])    
    %}
    %v=[mean(NC) max(NC) mean(edge_count) max(edge_count)];
    IData(i,:)=GraphImportance(G,i,Labels{i});
end

%% Compare COngestions Steiner vs Leaf

for i=2:3
    G=Gs{i};
    
    AP=find(or(G.Nodes.Type==1,G.Nodes.Type==3));
    edges_used=[];

    for s=1:length(AP)
        for t=s+1:length(AP)
           [P,d,eu] = shortestpath(G,AP(s),AP(t)); 
           edges_used=[edges_used eu];
        end
    end

    for e=1:size(G.Edges,1)
        edge_count(e) = sum(edges_used==e);
    end

    figure(258)
    subplot(1,2,i-1)
    pe = plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'NodeLabel',{},'EdgeAlpha',1,'LineWidth',edge_count./10);
    pe.EdgeCData = edge_count;
    colormap(flip(winter,1));
    title('Edges Congestion')

    %sgtitle(name)

end


%% Importance plot
INames={'Mean Node Betweeness','Max. Node Betweeness','Mean Edge Congestion','Max. Edge Congestion'};
figure(24)
hold on
for II=1:4
    subplot(2,2,II)
    %bar(IData(:,II)')
    bar(diag(IData(:,II)'),'stacked')
    set(gca,'xticklabel',Labels);
    title(INames{II})
end

sv=reshape(pw_dist,[16*16,3]);

figure(17)
boxplot(sv,'Notch','on','Labels',Labels)
title('Statistics of distance ratio')


% 1. Population
% 2. Buffer Area (2km)
% 3. Total Length
% 4. Sum dis_from_center
% 5. Buffer Density (pop/area)
% 6. Efficiency(?): (pop/total_length)
% 7. Mean ratio Path/Euclidean_dist -> From Center
% 8. Mean ratio Path/Euclidean_dist -> From every point
% 9. Buff Area/Total Length
plot_order=[1 2 3;
            6 9 5];
plot_titles={'Population', 'Buffer Area [km^{2}]','Network Length [km]';'Population per km of Network [hab/km]','Buffer Area per km of Network [km]','Population per unit area [hab/km^{2}]'};
Labels={'FN','ST','LV'};

%figure(18)
for p=1:2
    figure(18+p)
    for ff=1:3
    subplot(1,3,ff); 
    d = Data(plot_order(p,ff),1:3);
    bar(categorical(Labels),diag(d),'stacked')
    title(plot_titles(p,ff))
    end
    
end



%% Auxiliar Functions

function [v]=GraphImportance(G,i,name)
AP=find(or(G.Nodes.Type==1,G.Nodes.Type==3));
edges_used=[];

for s=1:length(AP)
    for t=s+1:length(AP)
       [P,d,eu] = shortestpath(G,AP(s),AP(t)); 
       edges_used=[edges_used eu];
    end
end

for e=1:size(G.Edges,1)
    edge_count(e) = sum(edges_used==e);
end

wbc = centrality(G,'betweenness','Cost',G.Edges.Weight);
n = numnodes(G);

brown=[173,96,10]./255;
c=copper(3);

figure(i)
subplot(1,3,1)
pn = plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'NodeLabel',{},'MarkerSize',(wbc+2).^0.4);
NC = 2*wbc./((n-2)*(n-1));
pn.NodeCData = NC;
colormap(flip(autumn,1));
title('Betweenness Centrality Scores - Weighted')

figure(i)
subplot(1,3,2)
pe = plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'NodeLabel',{},'EdgeAlpha',1,'LineWidth',edge_count./10);
pe.EdgeCData = edge_count;
colormap(flip(winter,1));
title('Edges Congestion')


figure(256)
subplot(1,3,i)
pe = plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'NodeLabel',{},'EdgeAlpha',1,'LineWidth',edge_count./10);
pe.EdgeCData = edge_count;
colormap(flip(winter,1));
title('Edges Congestion')

sgtitle(name)

%c = categorical({'apples','pears','oranges'});
%prices = [1.23 0.99 2.3];
%bar(c,prices)

wbc = centrality(G,'pagerank');

figure(i)
subplot(1,3,3)
pn = plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'NodeLabel',{},'MarkerSize',(wbc+2).^0.4);
NC = wbc;
pn.NodeCData = NC;
colormap(flip(autumn,1));
title('Betweenness PageRank Scores')

v=[mean(NC) max(NC) mean(edge_count) max(edge_count)];
end
%% Statistics for Density Layer

%{
Analyzed field: Density

Count: 1723

Unique values: 1720

NULL (missing) values: 0

Minimum value: 0.0

Maximum value: 21603.04118501

Range: 21603.04118501

Sum: 2504743.856753676

Mean value: 1453.7108861019594

Median value: 1095.21677932

Standard deviation: 1393.8213823628462

Coefficient of Variation: 0.9588023283641334

Minority (rarest occurring value): 18.4048088

Majority (most frequently occurring value): 0.0

First quartile: 734.206777435

Third quartile: 1682.8912814

Interquartile Range (IQR): 948.684503965

%}