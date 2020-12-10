%function cluster(edges)
%% INPUT VARIABLES
NX=50; %Nodes on X direction
NY=70; %Nodes on Y direction
D=1; %Lattice distance (side of equilateral triangle)
top=0.1;
%% IMPORT OUTPUTS

cases={'1','2','3','4'};

[~,FolderRoute] = uigetfile('*.txt','Chose Any File inside folder of domain');  % Opens Dialog to select txt that contains Nodes
DomainName='VSDN';

ncases=4;

MEI=cell(ncases,1);
MEF=cell(ncases,1);
Topo=cell(ncases,1);

for Case=1:ncases
    
        MEIstring=char(strcat('MEI-Case',char(cases(Case)),'-',DomainName,'.txt'));
        MEFstring=char(strcat('MEF-Case',char(cases(Case)),'-',DomainName,'.txt'));
        Topostring=char(strcat('Topo-Case',char(cases(Case)),'-',DomainName,'.txt'));
        
        MEI{Case}=dlmread(strcat(FolderRoute,MEIstring));
        MEF{Case}=dlmread(strcat(FolderRoute,MEFstring));
        Topo{Case}=dlmread(strcat(FolderRoute,Topostring));
        
        
end

dist={'A','B','C','D'};

nd=4;
Nodes=cell(nd,1);

for D=1:nd   
    NodeString=char(strcat(char(dist(D)),'-Nodes-',DomainName,'.txt'));        
    Nodes{D}=dlmread(strcat(FolderRoute,NodeString));
    [PP{D} CN{D}]=clusterFunction(MEF{D},Nodes{D},NX,NY,D,top);
end
%A) increase vs reo
figure(10)
hold on
for f=1:ncases
    scatter(MEI{f}(:,3),Topo{f}(:,3).*100)
end
legend('Case A','Case B','Case C')

%A) change of distribution

for Case=1:ncases
figure(10+Case)
hold on
hI = histogram(MEI{Case}(:,3));
dF=MEI{Case}(:,3).*(1+Topo{Case}(:,3));
%hF = histogram(MEF{Case}(:,3));
hF = histogram(dF);
%hC.Normalization = 'pdf';
hI.BinWidth = 0.025;
%hC.EdgeAlpha=0.1;
%hC.FaceAlpha=0.5;
%hT.Normalization = 'pdf';
hF.BinWidth = 0.025;
%hT.EdgeAlpha=0.1;
%hT.FaceAlpha=0.5;
legend('Initial','Final')
hold off
end

%Table
tableAT=zeros(3,4);
tableAPP=zeros(3,4);
midY=0.5*max(Topo{Case}(:,2));
for Case=1:ncases
    drT=mean(Topo{Case}(:,3));
    drU=mean(Topo{Case}((Topo{Case}(:,2)<midY),3));
    drD=mean(Topo{Case}((Topo{Case}(:,2)>=midY),3));
    tableAT(:,Case)=100.*[drT;drU;drD];
    
    drT=mean(Topo{Case}(PP{Case},3));
    CEUP=[];
    CEDN=[];
    for ppi=1:length(PP{Case})
        if Topo{Case}(PP{Case}(ppi),2)<midY
            CEUP(length(CEUP)+1)=Topo{Case}(PP{Case}(ppi),3);
        else
            CEDN(length(CEDN)+1)=Topo{Case}(PP{Case}(ppi),3);
        end
    end
    %indexUpper=Topo{Case}(PP{Case},2)<midY;
    %indexLower=Topo{Case}(PP{Case},2)>=midY;
    %drU=mean(Topo{Case}(indexUpper,3));
    %drD=mean(Topo{Case}(indexLower,3));
    drU=mean(CEUP);
    drD=mean(CEDN);
    tableAPP(:,Case)=100.*[drT;drU;drD];
    
end

%Pore volume Indexes from preferred path

TableB=zeros(3,3,ncases);
kv=4*pi()/3;

for Case=1:ncases    
    
    %TOTAL Domain
    FC=sum(sum(kv.*Nodes{Case}.^3)); %Total volumen of all the domain
    FU=sum(sum(kv.*Nodes{Case}(1:NY/2,:).^3)); %Total volumen of upper half of the domain
    FD=sum(sum(kv.*Nodes{Case}(1+NY/2:end,:).^3)); %Total volumen of lower (down) half of the domain
    
    %Connected domain
    
    CVT=0;
    CVU=0;
    CVD=0;
    for N=1:length(CN{Case}(:,1))
        if CN{Case}(N,1)<=(NY/2) %Upper domain
            CVT=CVT+ kv.*Nodes{Case}(CN{Case}(N,1),CN{Case}(N,2)).^3;
            CVU=CVU+ kv.*Nodes{Case}(CN{Case}(N,1),CN{Case}(N,2)).^3;
        else %Lower domain
            CVT=CVT+ kv.*Nodes{Case}(CN{Case}(N,1),CN{Case}(N,2)).^3;
            CVD=CVD+ kv.*Nodes{Case}(CN{Case}(N,1),CN{Case}(N,2)).^3;
        end
    end
    
    cp=100/FC;
    TableB(1,1,Case)=FC*cp;
    TableB(2,1,Case)=FU*cp;
    TableB(3,1,Case)=FD*cp;
    
    TableB(1,2,Case)=CVT*cp;
    TableB(2,2,Case)=CVU*cp;
    TableB(3,2,Case)=CVD*cp;
    
    TableB(1,3,Case)=CVT*100/FC;
    TableB(2,3,Case)=CVU*100/FU;
    TableB(3,3,Case)=CVD*100/FD;
    
end















