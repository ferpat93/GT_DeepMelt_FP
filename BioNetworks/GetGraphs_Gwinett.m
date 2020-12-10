% Gwinett Graphs

path='/Users/lfp3/Dropbox (GaTech)/GT/Spring-19/RootsPaper/ATL_map_files/Gwinett_WT_Dens_Centroids_N5.csv';
%path  = 'C:\Users\lfp3\Dropbox (GaTech)\GT\Spring-19\RootsPaper\ATL_map_files\Gwinett_WT_Dens_Centroids_N5.csv';
Points=dlmread(path,',',1,0);
a=Points(:,1:2); % a has to have just coordinates
points=a;

%{
GS=GetSteiner3(points);
GF=GetFan(points);
GL=GetLeaf(points);
%}

% Fan Network
gg=load('Gwinnett_FanGraph_N5.mat');
GF=gg.GF;
% Steiner Network:
gg=load('Gwinnett_SteinerGraph_N5.mat');
GS=gg.GS;
% Leaf Network
gg=load('Gwinnett_LeafGraph_N5.mat');
GL=gg.GL;

Gs={GF,GS,GL};
Labels={'Fan','Steiner','Leaf'};

write_lines=1;

if write_lines==1   
    for i=1:numel(Gs)

        path=fullfile(pwd,'ATL_map_files',strcat('Gwinnett_N5_Graph_Line_Layer_',Labels{i},'.csv'));

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











%% Auxiliar functions
% Fan
function [G]=GetFan(points)

nPoints=length(points(:,1));

Edges=ones(nPoints-1,2);
Edges(:,2)=transpose(linspace(2,nPoints,nPoints-1));
GraphNodes=[points 3*ones(nPoints,1)];
GraphNodes(1,3)=1;

% Build Graph

    % Assign Edge Lengths, width
Xc=GraphNodes(:,1);
Yc=GraphNodes(:,2);
xe=Xc(Edges);
ye=Yc(Edges);
weights=(sum([(xe(:,1)-xe(:,2)).^2 (ye(:,1)-ye(:,2)).^2],2)).^(1/2);

EdgeTable = table(Edges(:,1:2),weights,'VariableNames',{'EndNodes','Weight'});
NodeTable = table(GraphNodes(:,1),GraphNodes(:,2),GraphNodes(:,3),'VariableNames',{'X' 'Y','Type'});
G = graph(EdgeTable,NodeTable);

end

% Leaf
function [G]=GetLeaf(points)
scale=0.1; % OJO

points=points.*scale; % Scale points, largest scale: better resolution but exponentially increasing running time
nPoints=length(points(:,1));

%Simulation Parameters - 
%--> Set corners of the rectangular domain
xi=min(points(:,1));
xf=max(points(:,1));
yi=min(points(:,2));
yf=max(points(:,2));

%--> Set coordinates of petiole of the leaf
px=points(1,1);
py=points(1,2);

%--> Set growth rate (length of growth per iteration) and delta t
GR=1;
dt=1;

%Load Auxin Points
AuxPoints=[points(2:end,:) zeros(nPoints-1,1) ones(nPoints-1,1)];

%Set Kill Distance (multiples of GR - integers)
%% OJO -What does it do?
kd=2*GR; % Used to end the growth towards certain point 

%Creates struct object with indexed model
model=setModel(xi,xf,yi,yf,px,py,GR,dt,AuxPoints);
%display(model)

%Initializes nodes and edges
nodes=model.n0;
edges=[];
ConquerNodes=zeros(numel(model.AuxPoints(:,5)),1);

%Iteration over time steps - Until all Aux points have been reached
t=1;
while and(any(model.AuxPoints(:,5)==1),t<100000) % Until every point has been reached
    
    NewNodes=[]; % Empty matrix which will contain nodes created during the given time step
    [ActAP,model,ConquerNodes]=ActiveAP(t, model,nodes,kd,ConquerNodes); % Returns list of Active AP indexes based on time and space (time initially) 
    
    VorNodes=[nodes;model.nExt]; %Adds the synthetic external points for the Voronoi diagrams
    [VV,iV] = voronoin(VorNodes); % Create Voronin Diagrams for Vein Nodes each time

    AuxCoord=model.AuxPoints(ActAP,1:4);
        
    for VN=1:length(nodes(:,1)) % Loop along each node (excluding external)

       PolyV=VV(iV{VN},:); %Vertexes of polygon VN
       in = inpolygon(AuxCoord(:,1),AuxCoord(:,2),PolyV(:,1),PolyV(:,2)); % Aux points inside (in) each Voronoi diagram
       AttractingPoints=AuxCoord(in,:); %Asign points that attract each vein node (X and Y Coordinates)

       %Calculate new Node-Edge
       if isempty(AttractingPoints)==0 % if the array is not empty a new vein node is created
           mi=normr([ (AttractingPoints(:,1)-nodes(VN,1)) (AttractingPoints(:,2)-nodes(VN,2))]); % Finds the unitary vectors from the vein Node to each attracting Auxin Point
           mi=[mi(:,1).*AttractingPoints(:,4) mi(:,2).*AttractingPoints(:,4)]; % Multiplies each component by the weight of each auxin point
           nNodes=normr([sum(mi(:,1)) sum(mi(:,2))]); % Finds the norm of the resulting direction from the sum of all aux points (weighted)
           np=newvein(nNodes,nodes(VN,:));
           
           if any(ismember([nodes;NewNodes],np,'rows'))==0 %non-repeated node
            NewNodes=[NewNodes; np];
            edges=[edges;VN length(nodes(:,1))+length(NewNodes(:,1)) 0];
           end
           
       end 
       
    end

    nodes=[nodes;NewNodes]; %Adds the new nodes
    t=t+dt;
end

[edges,TrueEdges,TrueNodes]=veins(nodes,edges);

assignin('base','TrueEdges',TrueEdges);

grid=zeros(model.nx,model.ny);

for i=1:length(nodes(:,1))
   grid(nodes(i,1),nodes(i,2))=1; 
end

assignin('base','grid',grid);
sk= bwmorph(grid,'skel',Inf);

TrueNodes=[nodes(TrueNodes(:,1),1) nodes(TrueNodes(:,1),2) TrueNodes(:,2)];

%Gnodes=[model.n0 1]; %Adds stem as node
Gnodes = [];

[i,j] = find(bwmorph(sk,'branchpoints')); %Find Location of BranchPoints
Gnodes=[Gnodes; i j  2.*ones(length(i),1)]; %Adds Branching nodes

[i,j] = find(bwmorph(sk,'endpoints')); %Find Location of endPoints
Gnodes=[Gnodes; i j  3.*ones(length(i),1)]; %Adds EndPoints nodes


%DiffNodes=setdiff(TrueNodes(:,1:2),Gnodes(:,1:2),'rows'); %Compare two different approaches Gnodes usually simplifies


GraphNodes=Gnodes; % TrueNodes or Gnodes

for i=1:length(TrueEdges(:,1))  

    Ind_i=find(ismember(GraphNodes(:,1:2),nodes(TrueEdges(i,1),1:2),'rows'));
    Ind_f=find(ismember(GraphNodes(:,1:2),nodes(TrueEdges(i,2),1:2),'rows'));
    TrueEdges(i,1:2)=[Ind_i Ind_f] ; %Transform to GraphNodes indexes
    TrueEdges(i,4)= min(min(bwdistgeodesic(sk,GraphNodes(Ind_i,2),GraphNodes(Ind_i,1),'quasi')+ bwdistgeodesic(sk,GraphNodes(Ind_f,2),GraphNodes(Ind_f,1),'quasi')));
    if TrueEdges(i,4)==inf
        TrueEdges(i,4)=1;
    end
end





% Assign Edge Lengths, width
Edges=TrueEdges;
GraphNodes=[(xi+(GraphNodes(:,1)-1))./scale (yi+(GraphNodes(:,2)-1))./scale GraphNodes(:,3)];

Xc=GraphNodes(:,1);
Yc=GraphNodes(:,2);
xe=Xc(Edges(:,1:2));
ye=Yc(Edges(:,1:2));
weights=(sum([(xe(:,1)-xe(:,2)).^2 (ye(:,1)-ye(:,2)).^2],2)).^(1/2);

EdgeTable = table(TrueEdges(:,1:2),TrueEdges(:,3),TrueEdges(:,4), weights, ...
    'VariableNames',{'EndNodes' 'Width' 'Length' 'Weight'});
NodeTable = table(GraphNodes(:,1),GraphNodes(:,2),GraphNodes(:,3),'VariableNames',{'X' 'Y','Type'});
G = graph(EdgeTable,NodeTable);

brown=[173,96,10]./255;
c=copper(3);

figure(11)
plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'EdgeColor',brown,'NodeColor',c(G.Nodes.Type,:),'NodeLabel',{});
end
function [model]=setModel(xi,xf,yi,yf,px,py,GR,dt,AuxP)
%Creates Struct of mapped model based on the coordinates provided

%Number of X and Y positions (indexes)
model.nx=round(abs(xf-xi)/GR)+1;
model.ny=round(abs(yf-yi)/GR)+1;

%Index location of initial Vein Node
n0x=round(abs(px-xi)/GR)+1;
n0y=round(abs(py-yi)/GR)+1;
model.n0=[ n0x , n0y ];

nd=2; % Times of distance of the border

% External Nodes for the Voronoi polygons
model.nExt=[1-nd*(n0x) n0y;
            model.nx+nd*(model.nx-n0x) n0y;
            n0x 1-nd*(n0y);
            n0x model.ny+nd*(model.ny-n0y)];       
        
model.Xlim=[xi xf];
model.Ylim=[yi, yf];

AuxModel(:,1)=round(abs(AuxP(:,1)-xi)./GR)+1; % Bin X
AuxModel(:,2)=round(abs(AuxP(:,2)-yi)./GR)+1; % Bin Y
AuxModel(:,3)=round(AuxP(:,3)/dt)+1; % Time of appearance
AuxModel(:,4)=AuxP(:,4); % ??
AuxModel(:,5)=ones(length(AuxP(:,2)),1); % Boolean: If reached == 1

model.AuxPoints=AuxModel;
model.GR=GR;

end
function [ActAP,model,ConquerNode]=ActiveAP(t, model,nodes,kd,ConquerNode)

AuxPoints=model.AuxPoints;
ActAP=find(AuxPoints(:,5)==1); % Fifth column==1 means not reached yet

for i=1:length(ActAP)       
    [D,Node]=min(abs(AuxPoints(ActAP(i),1)-nodes(:,1))+abs(AuxPoints(ActAP(i),2)-nodes(:,2)));
    %if any(and(abs(AuxPoints(ActAP(i),1)-nodes(:,1))<kd,abs(AuxPoints(ActAP(i),2)-nodes(:,2))<kd))
    if D<=kd
        AuxPoints(ActAP(i),5)=0;   
        ConquerNode(ActAP(i))=Node;
        %disp(['Point reached:' num2str(ActAP(i)) ' at time: ' num2str(t)])
        %figure(10)
        %hold on
        %scatter(AuxPoints(ActAP(i),1),AuxPoints(ActAP(i),2))
      
    end
end

ActAP=intersect(find(AuxPoints(:,5)==1),find(AuxPoints(:,3)<=t)); % Finds indexes of not-absorbed -time-active AuxPoints
model.AuxPoints=AuxPoints;
end
function [nv]=newvein(v,c)
    x=v(1);
    y=v(2);
    t=y/x;
    
if t<=tand(25.5+90) % 3 or 7
    dx=0;
    if y>0 % 3
        dy=1;
    else % 7
        dy=-1;
    end     
elseif t<=tand(25.5+135) % 4 or 8
    if y>0 % 4
        dy=1;
        dx=-1;
    else % 8
        dy=-1;
        dx=1;
    end
elseif t<=tand(25.5) % 1 or 5
    dy=0;
    if x>0 % 1
        dx=1;
    else % 5
        dx=-1;
    end     
elseif t<=tand(25.5+45)  % 2 or 6
    if y>0 % 2
        dy=1;
        dx=1;
    else % 6
        dy=-1;
        dx=-1;
    end    
else
    dx=0;
    if y>0 % 3
        dy=1;
    else % 7
        dy=-1;
    end 
end

nv=[c(1)+dx c(2)+dy];

end
function [edges, TrueEdges, TrueNodes]=veins(nodes,edges)

a=1.5; %width assignment exponent

TrueEdges=[];

%Assign width to external edges
terminalNodes=setdiff(linspace(1,length(nodes(:,1)),length(nodes(:,1))),edges(:,1));
TrueNodes=[terminalNodes' 3.*ones(length(terminalNodes),1)];

AssignEdges=[];
for i=1:length(terminalNodes)
    idEdge=find(edges(:,2)==terminalNodes(i));
    edges(idEdge,3)=1; %assign initial width of 1 to terminal edges
    
    AssignEdges=[AssignEdges;[idEdge terminalNodes(i) edges(idEdge,3)]]; %Stores terminal edges as already assigned edges
    % Format [edgeID TrueEdge_EndNode Width]

end

while any(edges(:,3)==0) %Runs until every edge has a width
    
    New_AssignEdges=[];
    
    for i=1:length(AssignEdges(:,1))
        
        Parent_node=edges(AssignEdges(i,1),1);
        childs_i=find(edges(:,1)==Parent_node); %Edges that have the given node as father
        
        if length(childs_i)>1 %Is a branching Node (true)
            TrueNodes=unique([TrueNodes;[Parent_node 2]],'rows');
            TrueEdges=[TrueEdges;[Parent_node AssignEdges(i,2) AssignEdges(i,3)]];
            
            if all(edges(childs_i,3)>0) % All the widths of the childs are know, can calculate father's
                Parent_Edge=find(edges(:,2)==Parent_node); %Find the parent edge of all the childs
                
                if isempty(Parent_Edge)==0
                    edges(Parent_Edge,3)=(sum((edges(childs_i,3)).^a))^(1/a); %assigns its width
                    New_AssignEdges=[New_AssignEdges;[Parent_Edge Parent_node edges(Parent_Edge,3)]]; %Goes down one level
                end
            
            else
                New_AssignEdges=[New_AssignEdges;AssignEdges(i,:)];%Keeps edge on queue
            end
            
        else % Is an intermediate node
            Parent_Edge=find(edges(:,2)==Parent_node); %Find the parent edge of all the childs
            edges(Parent_Edge,3)=edges(childs_i,3); %assigns its width
            New_AssignEdges=[New_AssignEdges;[Parent_Edge AssignEdges(i,2) AssignEdges(i,3)]]; %Goes down one level            
        end              
    end
    AssignEdges=unique(New_AssignEdges,'rows');
    
end

for i=1:length(AssignEdges(:,1))
    TrueEdges=[TrueEdges;[1 AssignEdges(i,2) AssignEdges(i,3)]];
end

if any(TrueNodes(:,1)==1)
    Ind_1=find(TrueNodes(:,1)==1);
    TrueNodes(Ind_1,:)=[1 1];
    TrueNodes=unique(TrueNodes,'rows');
else
    TrueNodes=unique([TrueNodes;[1 1]],'rows');
end

TrueEdges=unique(TrueEdges,'rows');
end
    
% Steiner
function [G]=GetSteiner3(points)
%/Users/lfp3/Dropbox\ \(GaTech\)/GT/Previous\ Semesters/Summer-17/BioNetworksPaper/MPC15-master/experiments/stpToSmith.sh input.stp | /Users/lfp3/Dropbox\ \(GaTech\)/GT/Previous\ Semesters/Summer-17/BioNetworksPaper/MPC15-master/code/WarrenSmith/wds_smt_timing

nPoints=length(points(:,1));

%% Create .stp File
stpName='input.stp';
%outName='output.txt';

fout = fopen(stpName,'w');
strings=string({'SECTION Graph',['Nodes ',num2str(nPoints)],'END','SECTION Coordinates','END'});

% Start modifying text
for i=1:numel(strings)-1
    s=strings(i);
    fprintf(fout,'%s\n',s);
end

for i=0:nPoints-1
   s=['DD ', num2str(i) , ' ' , num2str(points(i+1,1)), ' ' , num2str(points(i+1,2))];
   fprintf(fout,'%s\n',s); 
end

fprintf(fout,'%s\n',strings{end}); 
fclose(fout);

%% Run Steiner Code
SCpath_terminal='/Users/lfp3/Dropbox\ \(GaTech\)/GT/Previous\ Semesters/Summer-17/BioNetworksPaper/MPC15-master/experiments';
SCpath_matlab='/Users/lfp3/Dropbox (GaTech)/GT/Previous Semesters/Summer-17/BioNetworksPaper/MPC15-master/experiments';
copyfile(stpName,SCpath_matlab)
%runString=['/SteinerExact -method Branch -input ', stpName];
%./stpToSmith.sh input.stp | ../code/WarrenSmith/wds_smt_timing
runString='./stpToSmith.sh input.stp | ../code/WarrenSmith/wds_smt_timing';
full_runString=['cd ' SCpath_terminal ' && ' runString]
[~,Out]=system(full_runString);
assignin('base','data',Out)

% Parse Output
[SteinerPoints, Edges]=parseOut2(Out);

% Remove replicate nodes

% Build Graph

GraphNodes=[points, 3.*ones(nPoints,1); SteinerPoints, 2.*ones(length(SteinerPoints(:,1)),1)];
GraphNodes(1,3)=1;

%Edges=Edges+1;

    % Assign Edge Lengths, width
Xc=GraphNodes(:,1);
Yc=GraphNodes(:,2);
xe=Xc(Edges);
ye=Yc(Edges);
weights=(sum([(xe(:,1)-xe(:,2)).^2 (ye(:,1)-ye(:,2)).^2],2)).^(1/2);

[SteinerPoints, Edges]=removeDuplicates(SteinerPoints,points,Edges,weights);


EdgeTable = table(Edges(:,1:2),weights,'VariableNames',{'EndNodes','Weight'});
    
NodeTable = table(GraphNodes(:,1),GraphNodes(:,2),GraphNodes(:,3),'VariableNames',{'X' 'Y','Type'});
G = graph(EdgeTable,NodeTable);


end
function [SteinerPoints, Edges]=parseOut(Out)

SP_i = strfind(Out,'point[');
SP_f = strfind(Out,'Edges:')-2;

Ed_i = SP_f+9;
Ed_f = strfind(Out,'Evaluated_Level:')-2;


SP_Char=strsplit(Out(SP_i:SP_f));
SteinerPoints=zeros(round(numel(SP_Char)/4),2);
for i=1:length(SteinerPoints(:,1))
    SteinerPoints(i,1)=str2double(SP_Char{i*4-1});
    SteinerPoints(i,2)=str2double(SP_Char{i*4});
end

Ed_Char=strsplit(Out(Ed_i:Ed_f));
Edges=zeros(round(numel(Ed_Char)/2),2);
for i=1:length(Edges(:,1))
    Edges(i,:)=[str2double(Ed_Char{i*2-1}) , str2double(Ed_Char{i*2})];
end

end
function [SteinerPoints, Edges]=parseOut2(Out)
Out = splitlines(Out);

Ee=find(contains(Out,'edges'));
El=Ee(end)+1;
SP=find(contains(Out,'steiner point coords'));
SP_i=SP(end)+1;
SP_f = El-2;
C = cellfun(@strsplit,Out(SP_i:SP_f),'UniformOutput',false);
SteinerPoints= [cell2mat(cellfun(@(x) str2double(x{2}),C,'UniformOutput',false)) cell2mat(cellfun(@(x) str2double(x{3}),C,'UniformOutput',false))];

Ed_Char=cellfun(@(x) strsplit(x,'-'),strsplit(Out{El},';'),'UniformOutput',false);
Edges = [cell2mat(cellfun(@(x) str2double(x{1}),Ed_Char(1:end-1),'UniformOutput',false)); cell2mat(cellfun(@(x) str2double(x{2}),Ed_Char(1:end-1),'UniformOutput',false))]';

end
function [SteinerPoints, Edges]=removeDuplicates(SteinerPoints,points,Edges,weights)

tol=1e-4*max(max(points(:,1))-min(points(:,1)),max(points(:,2))-min(points(:,2)));
fake_edges=Edges(find(weights<tol),:);
list=(1:(size(SteinerPoints,1)+size(points,1)))';
list(unique(fake_edges(:,2)))=[];
list=[(1:length(list))' list];
list(list(:,1)-list(:,2)==0,:)=[];

for i=1:size(fake_edges,1)
    Edges(Edges==fake_edges(i,2))=fake_edges(i,1);
end
Edges=unique(Edges,'rows');
Edges(Edges(:,1)-Edges(:,2)==0,:)=[];
SteinerPoints(fake_edges(:,2)-size(points,1),:)=[];

for i=1:size(list,1)
    Edges(Edges==list(i,2))=list(i,1);
end


end


