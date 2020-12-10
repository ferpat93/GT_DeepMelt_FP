%CREATES A-B-C-D ESCENARIOS FOR Nodes and Edges

function []=alternativecases(edges,nodes,mp,sp,mr,sr,fileName,FolderPath,NX,NY,D,kf)

format shortEng

k=2; %Scale constant
kn=k;
ke=3;
cnodes=coordNode(nodes,NX,NY,D);

%% Case A
% Edges
EdgesA=zeros(length(edges(:,1)),3); % Matrix for mean edge coordinates

for j=1:length(edges(:,1)) % Runs on every edge
    [x1,y1,x2,y2]=coord(j,NX,NY,D); %Gets initial and final point coordinates
    EdgesA(j,:)=[(x1+x2)*0.5 (y1+y2)*0.5 edges(j,1)];
end
EdgesRouteA=fullfile(FolderPath,strcat('HME-A-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case A
dlmwrite(EdgesRouteA,EdgesA);%writes named matrix to same folder

% Nodes
NodesRouteA=fullfile(FolderPath,strcat('HMN-A-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case A
dlmwrite(NodesRouteA,cnodes);%writes named matrix to same folder

%% Case D
% Edges
EdgesD=EdgesA; % Matrix for mean edge coordinates
EdgesflowD=edges;

for j=1:length(edges(:,1)) % Runs on every edge
    p= normcdf(edges(j,1),mr,sr); 
    newW=norminv(p,ke*mr,sr);
    EdgesD(j,3)=newW;
end

EdgesflowD(:,1)=EdgesD(:,3);
EdgesflowD(:,2)=kf*EdgesflowD(:,1).^4/D; %Sets the value alpha to the second row of the matrix   

EdgesRouteD=fullfile(FolderPath,strcat('D-Edges-',fileName)); % Sets name of Edges writen filename  
dlmwrite(EdgesRouteD,EdgesflowD);%writes named matrix to same folder

EdgesRouteD=fullfile(FolderPath,strcat('HME-D-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case D
dlmwrite(EdgesRouteD,EdgesD);%writes named matrix to same folder

% Nodes
cnodesD=cnodes; %Copies previous matrix

for j=1:length(cnodes(:,1)) % Runs on every edge
    p= normcdf(cnodes(j,3),mp,sp); 
    newW=norminv(p,kn*mp,sp);
    cnodesD(j,3)=newW;
end
NodesRouteD=fullfile(FolderPath,strcat('HMN-D-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case D
dlmwrite(NodesRouteD,cnodesD);%writes named matrix to same folder

NodesflowD(:,:,1)=nodes(:,:,1).*kn;

NodesRouteD=fullfile(FolderPath,strcat('D-Nodes-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case D
dlmwrite(NodesRouteD,NodesflowD);%writes named matrix to same folder
%% Case B and C

% Nodes
cnodesB=cnodes;
cnodesC=cnodes;

halfN=round(0.5*NY)*NX;

cnodesB(halfN+1:length(cnodes(:,1)),:)=cnodesD(halfN+1:length(cnodes(:,1)),:);
cnodesC(1:halfN,:)=cnodesD(1:halfN,:);

NodesRouteC=fullfile(FolderPath,strcat('HMN-C-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case C
dlmwrite(NodesRouteC,cnodesC);%writes named matrix to same folder
NodesRouteB=fullfile(FolderPath,strcat('HMN-B-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case B
dlmwrite(NodesRouteB,cnodesB);%writes named matrix to same folder

NodesflowB=nodes;
NodesflowB(round(0.5*NY)+1:NY,:,1)=NodesflowD(round(0.5*NY)+1:NY,:,1);

NodesRouteB=fullfile(FolderPath,strcat('B-Nodes-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case D
dlmwrite(NodesRouteB,NodesflowB);%writes named matrix to same folder

NodesflowC=nodes;
NodesflowC(1:round(0.5*NY),:,1)=NodesflowD(1:round(0.5*NY),:,1);

NodesRouteC=fullfile(FolderPath,strcat('C-Nodes-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case D
dlmwrite(NodesRouteC,NodesflowC);%writes named matrix to same folder


% Edges
cH=cnodesB(halfN,2);
EdgesB=EdgesA;
EdgesC=EdgesA;
for j=1:length(edges(:,1)) % Runs on every edge
   if EdgesB(j,2)>cH  
    EdgesB(j,3)=EdgesD(j,3);
   else
    EdgesC(j,3)=EdgesD(j,3);  
   end
end

EdgesRouteB=fullfile(FolderPath,strcat('HME-B-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case B
dlmwrite(EdgesRouteB,EdgesB);%writes named matrix to same folder
EdgesRouteC=fullfile(FolderPath,strcat('HME-C-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case C
dlmwrite(EdgesRouteC,EdgesC);%writes named matrix to same folder


EdgesflowC(:,1)=EdgesC(:,1);
EdgesflowC(:,2)=kf*EdgesflowC(:,1).^4/D; %Sets the value alpha to the second row of the matrix   

EdgesRouteC=fullfile(FolderPath,strcat('C-Edges-',fileName)); % Sets name of Edges writen filename  
dlmwrite(EdgesRouteC,EdgesflowC);%writes named matrix to same folder


EdgesflowB(:,1)=EdgesB(:,1);
EdgesflowB(:,2)=kf*EdgesflowB(:,1).^4/D; %Sets the value alpha to the second row of the matrix   

EdgesRouteB=fullfile(FolderPath,strcat('B-Edges-',fileName)); % Sets name of Edges writen filename  
dlmwrite(EdgesRouteB,EdgesflowB);%writes named matrix to same folder

