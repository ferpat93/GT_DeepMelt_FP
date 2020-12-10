% Load Data
Data=dlmread(fullfile(pwd,'ATL_map_files','Data_Table_ATL.txt'),'\t',1,0);
%[Population Area Perimeter Density Centroid_X Centroid_Y]

Center_X = 741594.849484141;
Center_Y = 3737918.386316658;
Center=[Center_X Center_Y];

Num_Clusters=5;
Num_Iter=100;
Coord=Data(:,5:6); % ** Change to meters?
Weights={Data(:,1), Data(:,4)} ; % Population

%% Order and Group
%order=[1 2 5 4 3]';
%group=[1 1 2 2 2]';


%%

% Regular K means
[Cluster_ID,Cluster_Coord] = kmeans(Coord,Num_Clusters);
path=fullfile(pwd,'ATL_map_files',strcat('NW_Centroids_N',num2str(Num_Clusters),'.csv'));
writeCSV(path,Cluster_Coord,Cluster_ID,Weights{1})

% Initially run regular K-means to initialize clusters
% Weighted K means

WType={'Pop','Dens'};

for i=1:2 % Loop Weight as density or population
    [ClusterID_W, Cluster_Coord_W, cluster_population, cluster_energy, it_num ] =  Weight_K_Means(Num_Clusters, Num_Iter , Coord , Weights{i} , Cluster_Coord );
    path=fullfile(pwd,'ATL_map_files',strcat('W_',WType{i},'_Centroids_N',num2str(Num_Clusters),'.csv'));
    writeCSV(path,Cluster_Coord_W,ClusterID_W,Weights{i})
end

%{
 ATL = shaperead('C:\Users\lfp3\Downloads\Zoning_District_Atl\Zoning_District.shp');
%}

%% Auxiliar Functions
function []=writeCSV(path,Coord,ID,pop)
    nC=size(Coord,1);    
    Weight=zeros(nC,1);
    X=Coord(:,1);
    Y=Coord(:,2);
    
    for ci=1:nC
        Weight(ci)=sum(pop(ID==ci));
    end

    Center_X = 741594.849484141;
    Center_Y = 3737918.386316658;
    %Center=[Center_X Center_Y];
    
    X=[Center_X; X];
    Y=[Center_Y; Y];
    Weight=[0; Weight];
    
    %T = table(X,Y,Weight,order,group);
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