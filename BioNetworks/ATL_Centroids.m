% Load Data
Data=dlmread(fullfile(pwd,'ATL_map_files','Data_Table_ATL.txt'),'\t',1,0);
%[Population Area Perimeter Density Centroid_X Centroid_Y]

Num_Clusters=15;
Num_Iter=1e16;
Coord=Data(:,5:6); % ** Change to meters?
Weights={Data(:,1), Data(:,4)} ; % Population - Density

% Regular K means
[Cluster_ID,Cluster_Coord] = kmeans(Coord,Num_Clusters);
path=fullfile(pwd,'ATL_map_files',strcat('NW_Centroids_N',num2str(Num_Clusters),'.csv'));
writeCSV(Cluster_Coord,Cluster_ID,Weights{1},path)

% Initially run regular K-means to initialize clusters
% Weighted K means

WType={'Pop','Dens'};

for i=1:1 % Loop Weight as density or population
    [ClusterID_W, Cluster_Coord_W, cluster_population, cluster_energy, it_num ] =  Weight_K_Means(Num_Clusters, Num_Iter , Coord , Weights{i} , Cluster_Coord );
    path=fullfile(pwd,'ATL_map_files',strcat('W_',WType{i},'_Centroids_N',num2str(Num_Clusters),'.csv'));
    writeCSV(Cluster_Coord_W,ClusterID_W,Weights{i},path)
end

u = repelem(Coord,Weights{1},1);
[Cluster_ID_WW,Cluster_Coord_WW] = kmeans(u,Num_Clusters);

wd=round(Weights{2}.*1000);
v = repelem(Coord,wd,1);
[Cluster_ID_WD,Cluster_Coord_WD] = kmeans(v,Num_Clusters);

figure()
scatter(Cluster_Coord(:,1),Cluster_Coord(:,2))
hold on
scatter(Cluster_Coord_WD(:,1),Cluster_Coord_WD(:,2),'*')
scatter(Cluster_Coord_WW(:,1),Cluster_Coord_WW(:,2),'+')

%{
 ATL = shaperead('C:\Users\lfp3\Downloads\Zoning_District_Atl\Zoning_District.shp');
%}

%% Auxiliar Functions
function []=writeCSV(Coord,ID,pop,path)
    nC=size(Coord,1);    
    Weight=zeros(nC,1);
    X=Coord(:,1);
    Y=Coord(:,2);
    
    for ci=1:nC
        Weight(ci)=sum(pop(ID==ci));
    end

    T = table(X,Y,Weight);
    writetable(T,path,'WriteRowNames',true) 
end



%{

function  [d]=HaversineDist(lat,long)

    long=deg2rad(long);
    lat=deg2rad(lat);
    dlat=diff(lat);
    dlong=diff(long);
    a=sin(dlat/2).^2+cos(lat(1:end-1)).*cos(lat(2:end)).*sin(dlong/2).^2;
    c=2*atan2(sqrt(a),sqrt(1-a));
    R=6371000; %in metres
    d=R.*c;
  
end

%}