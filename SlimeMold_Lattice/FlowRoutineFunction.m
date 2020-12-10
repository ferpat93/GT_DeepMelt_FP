%% Flow Routine FUNCTION

function [fs]=FlowRoutineFunction(Path,NameEdges,NameNodes,caseString)

gplot=0; %Generate plot? true=1 false=0

m=strcat(Path,NameEdges);
edges=dlmread(m);     % Stores as 'nodes' the  matrix
edges0=edges;

%[fileNameT,FolderPathT] = uiputfile('*.txt','Save Final oputput as: ');  % Opens Dialog to select location of saved files
fileName=NameEdges;
fileName(1:8) = [];
fileName(end-3:end)=[];

filepath2=fullfile(Path,NameNodes); % Set imported file name
nodes2D=importdata(filepath2);      % Stores as 'nodes2D' the  matrix
        
NY=length(nodes2D(:,1,1)); %Computex NY from nodes matrix
NX=(length(edges(:,1))+2*NY-1)/(3*NY-2); %computes NX from total nodes number
NAN=length(nodes2D(1,:))/NX; %Number of atributes of nodes matrix (length of 3rd dimension)

nodes=reshape(nodes2D,[NY NX NAN]); %Turns Nodes2D into a 3D matrix 'nodes'
clear nodes2D % Erases value of nodes2D to save memory

%D = input('Type D value:  '); % asks the user for the lattice fundamental length D
D=1;

me=0.480;%e-6; %Mean of edge radii
de=0.120;%e-6; %Std deviation of edge radii
%% Sampling technique - is replacing the selection

%NS=ceil(0.3*NX*NY); % Number of nodes chosen to induce flow
NS=1;
VectorNodes=reshape(nodes(:,:,1),[NX*NY 1]); % Creates a vector containing the PoreSize of each Node - weight of the random sampling
sampleNodes = randsample(NX*NY,NS,true,VectorNodes); %WEIGHTED Creates a vector containing the positions on VectorNodes which will be used
%sampleNodes = randsample(NX*NY,NS,true); %NON WEIGHTED Creates a vector containing the positions on VectorNodes which will be used

sources=zeros(NS,2); % Matrix of sources to compute 
sources(:,2)=ceil(sampleNodes/NY); %Finds row index (Y direction)
sources(:,1)=sampleNodes-(sources(:,2)-1)*NY; %Finds column index (X direction)

if gplot==1
%Number of plots to be created
Nplot=2; %Number of plots that will be created
vplot=zeros(Nplot,1); %vector containing the source numbers where the plot will be built
for np=1:Nplot-1 
    vplot(np)=np*round(length(sources(:,1))/Nplot); %Sets the source where plot will be performed - evenly spaced along sources with Nplot steps
end
vplot(Nplot)=length(sources(:,1)); %Last plot after the last source
%[fileName,FolderPath] = uiputfile('*.fig','Location and name of saved figure');  % Opens Dialog to select location of saved files
end

%% Generation of flow

is=1; %Sets the indexes for the sink point (first row)
js=round(NX/2); %Sets the indexes for the sink point (mid axis x)
sink=(is-1)*NX+js; %ID of the sink based on row and column

FL=zeros(length(sources(:,1)),1);
%HA=zeros(length(sources(:,1)),5); %Vector of min max and percentiles of
%qEdge
%cellplot=cell(length(sources(:,1)),1); %CELL ARRAY FOR THE PLOT

for s=1:length(sources(:,1)) % Iterate over every source   
    
    [aP]=createAP(NX,NY,edges); % Creation of matrix of hydraulic pressure
    
    %% Generation of flow
    %Sets qin and qout values in the corresponding node
    q=zeros(NX*NY,1); % Solution vector containing q   
    i=sources(s,1); %Gets the indexes of evaluated sink
    j=sources(s,2); %Gets the indexes of evaluated sink
    sourceID=(i-1)*NX+j; %ID of the node based on row and column
    q(sourceID)=(-4*pi()/3)*nodes(i,j,1).^3; % Sets the given stream value for the given source
    q(sink)=-q(sourceID); %Sets the output stream at the sink
    %Solves the matrix to fet the vector 
    %P=(aP)\(q); %Solves the set of linear equations to find Pressure values for each node
    %[P,flag] = symmlq(aP,q,1e-3,50); %Solves the set of linear equations to find Pressure values for each node
    %FL(s)=flag;
    %[P] = pcg(aP,q,1e-11,NX*NY*3); %Solves the set of linear equations to find Pressure values for each node
    L = ichol(sparse(aP),struct('michol','on'));
    [P,flag] = pcg(aP,q,1e-8,NX*NY*3,L,L'); %Solves the set of linear equations to find Pressure values for each node
    FL(s)=flag;
    %Compute flow on each edge
    qEdge=zeros(length(edges(:,1)),1); % Vector containing flow value for each edge
    for nedge=1:length(edges(:,1)) % Iterates through each edge
        [i1,j1,i2,j2]=source(nedge,NX,NY); %Gets the id of the neighbor edge and then the indexes of the underlying  nodes
        ID1=(i1-1)*NX+j1; %ID of the node based on row and column
        ID2=(i2-1)*NX+j2; %ID of the node based on row and column
        qEdge(nedge,1)=edges(nedge,2)*(P(ID1)-P(ID2)); %Assigns the flow value to the node based on alpha(edge) and pressure difference between nodes     
    end
    
    edges=updateEdges2(D,edges,abs(qEdge),q(sourceID),me,de,NS); %Updates the mesh properties calling the auxiliar function
    
    %qEdge Analysis
    %histogram(abs(qEdge),50)
    %HA(s,:)=[min(abs(qEdge)) prctile(abs(qEdge),25) prctile(abs(qEdge),50) prctile(abs(qEdge),75) max(abs(qEdge))]./(0.01*q(sink));
    %cellplot{s}=[edges(:,2) abs(qEdge)./(0.01*q(sink))];
    
    if gplot==1 && ismember(s,vplot)==1
    spc=round(100*(s/length(sources(:,1))));
    FigRoute=fullfile(dirname,strcat(num2str(spc),'%-',fileName,'.fig')); % Sets name of Nodes writen filename 
    plotmesh(1,edges,nodes,NX,NY,D,FigRoute); %Calls the auxiliar function that plots and prints figure
    end 
    
end

fs=sum(FL);

topology=dW(edges./edges0,NX,NY,D);
MEI=dW(edges0,NX,NY,D);
MEF=dW(edges,NX,NY,D);
topoRoute=fullfile(Path,strcat('Topo-',caseString)); % Sets name of Edges writen filename  
MEFRoute=fullfile(Path,strcat('MEF-',caseString)); % Sets name of Edges writen filename  
MEIRoute=fullfile(Path,strcat('MEI-',caseString)); % Sets name of Edges writen filename  
dlmwrite(topoRoute,topology);%writes named matrix to same folder
dlmwrite(MEFRoute,MEF);%writes named matrix to same folder
dlmwrite(MEIRoute,MEI);%writes named matrix to same folder
