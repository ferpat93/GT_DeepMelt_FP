%% PLOT MESH
function []=plotmesh(save,edges,nodes,NX,NY,D,FigRoute)

LY=(NY-0.5)*D; % length of domain in X axis
LX=(NX-1)*cos(deg2rad(30))*D; %length of domain in Y axis
k=max(edges(:,1));

%figure(1) %Lattice mesh
figure
clf %resets the figure
ax=axes; %names axes of figure
set(ax,'Ydir','reverse') % Set Y axis on reverse direction
axis([0 LX 0 LY]) %Sets axis boundaries
axis equal % Sets the same scale on both axis
hold on % Holds the assigment of data sets to the same figure - each link is one data set
for j=1:length(edges(:,1)) % Runs on every edge
    [x1,y1,x2,y2]=coord(j,NX,NY,D); %Gets initial and final point coordinates
    L1 = plot([x1;x2],[y1;y2],'color',[.3 .4 0.6]); %Plots the edge and assigns color
    set(L1,'LineWidth',(4*edges(j,1)/k+0.001)) %Sets linewidth according to r value
end

cn=coordNode(nodes,NX,NY,D);%Generate matrix with nodes coordinates and pore size (oil)
scatter(cn(:,1),cn(:,2),10,cn(:,3),'filled'); %Plot Nodes with size on colormap

if save==1                                                               
  %  print(figure(1),FigRoute,'-dtiffn'); %Code for saving the figure
    saveas(figure(1),FigRoute,'fig'); %Code for saving the figure
end

