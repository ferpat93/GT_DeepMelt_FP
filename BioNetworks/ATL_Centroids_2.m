% Load Data
Data=dlmread(fullfile(pwd,'ATL_map_files','LayerFiles','Geo_data_meters.csv'),',',1,0);
%[Population Area Density Centroid_X Centroid_Y]

Num_Clusters=10;
Num_Iter=1e16;
Coord=Data(:,4:5); % ** Change to meters?
Weights={Data(:,1), Data(:,3)} ; % Population - Density

% Regular K means
[Cluster_ID,Cluster_Coord] = kmeans(Coord,Num_Clusters);
path=fullfile(pwd,'ATL_map_files',strcat('NW_Centroids_N',num2str(Num_Clusters),'.csv'));
writeCSV(Cluster_Coord,Cluster_ID,Weights{1},path)

% Initially run regular K-means to initialize clusters
% Weighted K means

WType={'Pop','Dens'};

%% Use Weighted K_Means
%{
for i=1:1 % Loop Weight as density or population
    [ClusterID_W, Cluster_Coord_W, cluster_population, cluster_energy, it_num ] =  Weight_K_Means(Num_Clusters, Num_Iter , Coord , Weights{i} , Cluster_Coord );
    path=fullfile(pwd,'ATL_map_files',strcat('W_',WType{i},'_Centroids_N',num2str(Num_Clusters),'.csv'));
    writeCSV(Cluster_Coord_W,ClusterID_W,Weights{i},path)
end
%}

%% Use tricked-traditional k means

% Population
W=Weights{1};
u = repelem(Coord,W,1);
[ClusterID_W, Cluster_Coord_W ] = kmeans(u,Num_Clusters); 
path=fullfile(pwd,'ATL_map_files',strcat('WT_',WType{1},'_Centroids_N',num2str(Num_Clusters),'.csv'));
d=unique([u ClusterID_W],'rows');
clearvars u ClusterID_W
writeCSV(Cluster_Coord_W,d(:,3),Weights{1},path)

% Density
W=round(Weights{2}.*10);
u = repelem(Coord,W,1);
[ClusterID_W, Cluster_Coord_W ] = kmeans(u,Num_Clusters);
path=fullfile(pwd,'ATL_map_files',strcat('WT_',WType{2},'_Centroids_N',num2str(Num_Clusters),'.csv'));
d=unique([u ClusterID_W],'rows');
clearvars u ClusterID_W
writeCSV(Cluster_Coord_W,d(:,3),Weights{1},path)

%{
figure()
scatter(Cluster_Coord(:,1),Cluster_Coord(:,2))
hold on
scatter(Cluster_Coord_WD(:,1),Cluster_Coord_WD(:,2),'*')
scatter(Cluster_Coord_WW(:,1),Cluster_Coord_WW(:,2),'+')
%}


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
    
    Center_X = 741594.849484141;
    Center_Y = 3737918.386316658;
    %Center=[Center_X Center_Y];
    
    X=[Center_X; X];
    Y=[Center_Y; Y];
    Weight=[0; Weight];
    
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