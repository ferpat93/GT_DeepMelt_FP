%CREATES A-B-C-D ESCENARIOS FOR Nodes and Edges

function []=alternativecasesVsdn(~,nodes,mp,sp,~,~,fileName,FolderPath,NX,NY,D,~)

format shortEng

k=2; %Scale constant
kn=k;
ke=1;
cnodes=coordNode(nodes,NX,NY,D);

%% Case A
% Nodes
NodesRouteA=fullfile(FolderPath,strcat('HMN-A-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case A
dlmwrite(NodesRouteA,cnodes);%writes named matrix to same folder

%% Case D
% Nodes
cnodesD=cnodes; %Copies previous matrix

for j=1:length(cnodes(:,1)) % Runs on every edge
    p= normcdf(cnodes(j,3),mp,sp); 
    newW=norminv(p,mp,mp);
    cnodesD(j,3)=newW;
end

NodesRouteD=fullfile(FolderPath,strcat('HMN-100-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case D
dlmwrite(NodesRouteD,cnodesD);%writes named matrix to same folder

NodesflowD(:,:,1)=nodes(:,:,1).*kn;

NodesRouteD=fullfile(FolderPath,strcat('100-Nodes-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case D
dlmwrite(NodesRouteD,NodesflowD);%writes named matrix to same folder
%% Case B and C

%% Case D
% Nodes
cnodesD=cnodes; %Copies previous matrix

for j=1:length(cnodes(:,1)) % Runs on every edge
    p= normcdf(cnodes(j,3),mp,sp); 
    newW=norminv(p,mp,mp*0.05);
    cnodesD(j,3)=newW;
end

NodesRouteD=fullfile(FolderPath,strcat('HMN-005-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case D
dlmwrite(NodesRouteD,cnodesD);%writes named matrix to same folder

NodesflowD(:,:,1)=nodes(:,:,1).*kn;

NodesRouteD=fullfile(FolderPath,strcat('005-Nodes-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case D
dlmwrite(NodesRouteD,NodesflowD);%writes named matrix to same folder

%% Case D
% Nodes
cnodesD=cnodes; %Copies previous matrix

for j=1:length(cnodes(:,1)) % Runs on every edge
    p= normcdf(cnodes(j,3),mp,sp); 
    newW=norminv(p,mp,mp*0.5);
    cnodesD(j,3)=newW;
end

NodesRouteD=fullfile(FolderPath,strcat('HMN-50-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case D
dlmwrite(NodesRouteD,cnodesD);%writes named matrix to same folder

NodesflowD(:,:,1)=nodes(:,:,1).*kn;

NodesRouteD=fullfile(FolderPath,strcat('50-Nodes-',fileName)); % Sets name of matrix for heatmap of oilconcentration -  case D
dlmwrite(NodesRouteD,NodesflowD);%writes named matrix to same folder
