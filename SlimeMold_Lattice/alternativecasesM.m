%CREATES A-B-C-D ESCENARIOS FOR Nodes and Edges

function []=alternativecasesM(edges,nodes,mp,sp,mr,sr,fileName,FolderPath,NX,NY,D,kf)

format shortEng

k=1.5; %Scale constant
kn=k;
ke=k;
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
%% Case B

% Edges
EdgesB=EdgesA; % Matrix for mean edge coordinates
EdgesflowB=edges;

for j=1:length(edges(:,1)) % Runs on every edge
    p= normcdf(edges(j,1),mr,sr); 
    newW=norminv(p,kmb*mr,ksb*sr);
    EdgesD(j,3)=newW;
end

EdgesflowB(:,1)=EdgesB(:,3);
EdgesflowB(:,2)=kf*EdgesflowB(:,1).^4/D; %Sets the value alpha to the second row of the matrix   

EdgesRouteB=fullfile(FolderPath,strcat('B-Edges-',fileName)); % Sets name of Edges writen filename  
dlmwrite(EdgesRouteB,EdgesflowB);%writes named matrix to same folder

EdgesRouteB=fullfile(FolderPath,strcat('HME-B-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case D
dlmwrite(EdgesRouteB,EdgesB);%writes named matrix to same folder

% Nodes
cnodesB=cnodes; %Copies previous matrix

for j=1:length(cnodes(:,1)) % Runs on every edge
    p= normcdf(cnodes(j,3),mp,sp); 
    newW=norminv(p,k*mp,sp);
    cnodesD(j,3)=newW;
end
NodesRouteD=fullfile(FolderPath,strcat('HMN-D-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case D
dlmwrite(NodesRouteD,cnodesD);%writes named matrix to same folder

NodesflowD(:,:,1)=nodes(:,:,1).*kn;

NodesRouteD=fullfile(FolderPath,strcat('D-Nodes-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case D
dlmwrite(NodesRouteD,NodesflowD);%writes named matrix to same folder
