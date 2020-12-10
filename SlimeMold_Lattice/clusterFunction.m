function [pp, CN]=clusterFunction(MEF,nodes,NX,NY,D,top)
%% INPUT VARIABLES
%NX=50; %Nodes on X direction
%NY=70; %Nodes on Y direction
%D=1; %Lattice distance (side of equilateral triangle)
%top=0.1;
%[pathname,FolderSource] = uigetfile('*.txt','Choose edges file');  % Opens Dialog to select txt that contains Nodes
%m=strcat(FolderSource,pathname);
%edges=dlmread(m);     % Stores as 'nodes' the  matrix

%[pathname,FolderSource] = uigetfile('*.txt','Choose nodes file');  % Opens Dialog to select txt that contains Nodes
%m=strcat(FolderSource,pathname);
%nodes=dlmread(m);     % Stores as 'nodes' the  matrix
%%
[~,I]=sort(MEF(:,3),'descend'); %sorts matrix of edges by highest to minimum radius
      
%edges(:,2)=edgesP(:,1);
%edges(:,1)=linspace(1,length(edges(:,2)),length(edges(:,2)));

remedges=I(1:round(length(I)*top)); %Remaining edges are top according to the criteria

idCluster=0;
Clusters={};

while isempty(remedges)==0 %Runs until every edge is assigned to a cluster
    
    idCluster=idCluster+1;  %Gets here when edges to connect is empty -  new cluster
    EdgestoConnect=remedges(1);
    remedges(1)=[];
    Clusters{idCluster}=[];
    
    while isempty(EdgestoConnect)==0 %Do until no more edges of the cluster need to be connected
        Clusters{idCluster}=[Clusters{idCluster} ; EdgestoConnect];
        [EdgestoConnect,remedges]=edgeNeighbours(EdgestoConnect,NX,NY,remedges);
    end

end

[vm, im]=max(cellfun(@length,Clusters));
vm
%% Index 3
% Volume of pores connected by preferred path

pp=Clusters{im};
[i1,j1,i2,j2]=source(pp(1),NX,NY);
CN=[[i1 j1 nodes(i1,j1)] ; [i2 j2 nodes(i2,j2)]];
for ne=2:length(pp)
[i1,j1,i2,j2]=source(pp(ne),NX,NY);
CN=[CN; [i1 j1 nodes(i1,j1)] ; [i2 j2 nodes(i2,j2)]];
end

CN=unique(CN,'rows'); %CN connected nodes [indexi indexj r]


%Show 
showhistogram=0;
if showhistogram==1
figure(10)
hold on
hC = histogram(CN(:,3));
hT = histogram(nodes);
hC.Normalization = 'pdf';
hC.BinWidth = 0.5;
%hC.EdgeAlpha=0.1;
%hC.FaceAlpha=0.5;
hT.Normalization = 'pdf';
hT.BinWidth = 0.5;
%hT.EdgeAlpha=0.1;
%hT.FaceAlpha=0.5;
legend('Full Network','Connected path')
hold off
end




%% PLOT Preferred path

save=0;
index=pp;

LY=(NY-0.5)*D; % length of domain in X axis
LX=(NX-1)*cos(deg2rad(30))*D; %length of domain in Y axis
k=max(MEF(:,3));

%figure(1) %Lattice mesh
figure
clf %resets the figure
ax=axes; %names axes of figure
set(ax,'Ydir','reverse') % Set Y axis on reverse direction
axis([0 LX 0 LY]) %Sets axis boundaries
%axis equal % Sets the same scale on both axis
hold on % Holds the assigment of data sets to the same figure - each link is one data set
for j=1:length(index) % Runs on every edge
    [x1,y1,x2,y2]=coord(index(j),NX,NY,D); %Gets initial and final point coordinates
    L1 = plot([x1;x2],[y1;y2],'color',[.3 .4 0.6]); %Plots the edge and assigns color
    set(L1,'LineWidth',(4*MEF(index(j),3)/k+0.001)) %Sets linewidth according to r value
end


if save==1                                                               
  %  print(figure(1),FigRoute,'-dtiffn'); %Code for saving the figure
    saveas(figure(1),FigRoute,'fig'); %Code for saving the figure
end
