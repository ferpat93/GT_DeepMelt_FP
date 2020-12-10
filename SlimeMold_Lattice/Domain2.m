%Set domain of lattice
%Two dimensions
global PEDGE PNODE

%% INPUT VARIABLES
NX=40; %Nodes on X direction
NY=50; %Nodes on Y direction
D=1; %Lattice distance (side of equilateral triangle)

mp=1.45;%e-14; %Mean of pore diameter %values in nanometers
%sp=mp*0.75;
sp=0.35;%e-14; %Std deviation of pore diameter

mr=0.480;%e-12; %Mean of edge radii
%sr=mr*0.000001;
sr=0.120;%e-12; %Std deviation of edge radii

n=1;%e-3; %Viscosity of the fluid
kf=pi()/(8*n); % Constant of flow

gplot=0; %Generate plot? true=1 false=0
%% SET EDGES

NE=3*NX*NY-2*(NX+NY)+1; % Number of Lattice edges from number of nodes

edges=zeros(NE,2); %Dimension to be decided - Matrix of dges attributes
NAN=1; %number of attributes to be stored on node matrix, length of 3rd dimension
nodes=zeros(NY,NX,NAN);% Matrix of nodes and its properties

%Sets edge width
for k=1:length(edges(:,1))
    edges(k,1)=norminv(PEDGE(k),mr,sr); %set normally distributed width of edge  
    if edges (k,1)<0 %if the width is negative assumes closed fracture
        edges(k,1)=0; %assigns closed value
    end    
end

edges(:,2)=kf*edges(:,1).^4/D; %Sets the value alpha to the second row of the matrix

%Sets pore size
for k=1:NX %Iterates over X
    for f=1:NY %Iterates over Y
    nodes(f,k,1)=norminv(PNODE(f,k),mp,sp); %set normally distributed size of pore
        if nodes(f,k,1)<0 %if the width is negative assumes closed fracture
            nodes(f,k,1)=0; %assigns closed value
        end
    end
end

[fileName,FolderPath] = uiputfile('*.txt','Save domain ');  % Opens Dialog to select location of saved files
NameDomain=fileName;
NameDomain(end-3:end)=[];
mkdir(FolderPath,NameDomain)
FolderPath=strcat(FolderPath,NameDomain);
NodesRoute=fullfile(FolderPath,strcat('A-Nodes-',fileName)); % Sets name of Nodes writen filename                                                                
EdgesRoute=fullfile(FolderPath,strcat('A-Edges-',fileName)); % Sets name of Edges writen filename                                                      
dlmwrite(EdgesRoute,edges); %writes in the destination folder the 2D matrix with node information

nodes2D=reshape(nodes,[NY NX*NAN]); %Creates 2D matrix from nodes, with NY rows and (NAN*NX) columns
dlmwrite(NodesRoute,nodes2D);%writes node2D matrix into a .txt file on the destination folder

%alternativecases(edges,nodes,mp,sp,mr,sr,fileName,FolderPath,NX,NY,D,kf)
%% PLOT MESH

if gplot==1 %if plotting and saving the mesh is needed
   %[fileName,FolderPath] = uiputfile('*.tiff','Save domain ');  % Opens Dialog to select location of saved files
   %[fileName,FolderPath] = uiputfile('*.fig','Set figure location ');  % Opens Dialog to select location of saved files
   FigRoute=fullfile(FolderPath,strcat(NameDomain,'-Initial.fig')); % Sets name of Nodes writen filename   
   plotmesh(1,edges,nodes,NX,NY,D,FigRoute) %Auxiliar function that plotes and prints the figure
end

%%

