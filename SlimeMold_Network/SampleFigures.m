close all
colormap jet
Nodes = [ 2 3 1 2 3 4 1 2 3 4 2 3; 0 0 2 2.0001 2.0001 2 5 5 5 5 7 7]';
Euc_Dist = triu(squareform(pdist(Nodes)));

Edges = [1 2; 1 4; 2 5; 3 4; 4 5; 5 6; 4 8; 5 9; 7 8; 8 9; 9 10; 8 11; 9 12; 11 12]; 
Width = ones(size(Edges,1),1);
Length = ones(size(Edges,1),1);
for i=1:size(Edges,1)
    Length(i) = Euc_Dist(Edges(i,1),Edges(i,2));
end
 
NodeTable = table(Nodes(:,1),Nodes(:,2),'VariableNames',{'X' 'Y'});
EdgeTable = table(Edges,Width,Length,'VariableNames',{'EndNodes' 'Width' 'Length'});
G = graph(EdgeTable,NodeTable);
GD = G;
Width(8) = 4; % Width(13) = 4;
GD.Edges.Width = Width;
G.Edges.drag = (G.Edges.Length)./(G.Edges.Width).^4;
GD.Edges.drag = (GD.Edges.Length)./(GD.Edges.Width).^4;
G = betweenness(G); GD = betweenness(GD);

GG = GD ; GG.Edges.Weight = GD.Edges.drag;
MST = betweenness(minspantree(GG));

tri = delaunay(Nodes);
A = adjacency(digraph(tri, tri(:, [2 3 1]))); A = A | A';
w = A.*(Euc_Dist+Euc_Dist');
DT = graph(w);
DT.Edges.Length = DT.Edges.Weight;
DT.Edges.Width = ones(size(DT.Edges,1),1);
DT.Edges.drag = (DT.Edges.Length)./(DT.Edges.Width).^4;
DT.Nodes.X = Nodes(:,1); DT.Nodes.Y = Nodes(:,2); 
DT = betweenness(DT);

% figure; hold on
% scatter(Nodes(:,1),Nodes(:,2))
% triplot(tri,Nodes(:,1),Nodes(:,2));

figure;
subplot(1,4,4);
% plot(DT,'XData',DT.Nodes.X,'YData',DT.Nodes.Y,'EdgeColor',[0 0 0.3],'NodeColor',[1 0 0]);
% xlabel('X Coordinate'); ylabel('Y Coordinate'); axis equal; axis tight
plotGraph(DT); title('DT')
subplot(1,4,3);
% plot(MST,'XData',MST.Nodes.X,'YData',MST.Nodes.Y,'EdgeColor',[0 0 0.3],'NodeColor',[1 0 0]);
% xlabel('X Coordinate'); ylabel('Y Coordinate'); axis equal; axis tight
plotGraph(MST); title('MST')
subplot(1,4,2); 
plotGraph(GD); title('G')
% plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'EdgeColor',[0 0 0.3],'NodeColor',[1 0 0],'LineWidth',G.Edges.Width);
% xlabel('X Coordinate'); ylabel('Y Coordinate'); axis equal; axis tight
subplot(1,4,1);
plotGraph(G); title('U')

colormap(jet);
caxis([0 0.65])
colorbar('southoutside')

%% Auxiliar

function []=plotGraph(G)
    %map = brewermap(101,'YlOrRd');
    map = jet(101);
    maxv = 0.65;
    edgePropToColor = G.Edges.BD;
    scaledV= round(100*edgePropToColor./maxv);
    scaledV(scaledV>100)=100;
    
    NodeC=brewermap(1,'Set1');
    Colors = map(scaledV+1,:);
    %figure;
    %plot(G,'EdgeColor',Colors,'NodeColor',NodeC(G.Nodes.Type+1,:),'NodeLabel',{});
    P = plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'EdgeColor',Colors,'NodeColor',NodeC,'LineWidth',G.Edges.Width,'NodeFontAngle','italic','NodeFontWeight','bold');
    P.EdgeAlpha = 1;
    xlabel('x')
    ylabel('y')
    axis equal
    axis tight

end

function G = betweenness(G)
    AP = 1:G.numnodes;
    GD = G; 
    G.Edges.Weight = G.Edges.Length; %PP = triu(distances(G));
    GD.Edges.Weight = G.Edges.drag; %PD = triu(distances(GD));
    
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

    G.Edges.BD = (2./((numel(AP)-2)*(numel(AP)-1))).*edge_count_D;
    G.Edges.BL = (2./((numel(AP)-2)*(numel(AP)-1))).*edge_count_P;

end