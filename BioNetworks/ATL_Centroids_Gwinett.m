% Load Data
Data=dlmread(fullfile(pwd,'ATL_map_files','LayerFiles','Geo_data_meters_Gwinett.csv'),',',1,0);
%[Population Area Density Centroid_X Centroid_Y]

Num_Clusters=15;
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


%% Use tricked-traditional k means

% Population
W=Weights{1};
u = repelem(Coord,W,1);
[ClusterID_W, Cluster_Coord_W ] = kmeans(u,Num_Clusters); 
path=fullfile(pwd,'ATL_map_files',strcat('Gwinett_WT_',WType{1},'_Centroids_N',num2str(Num_Clusters),'.csv'));
d=unique([u ClusterID_W],'rows');
clearvars u ClusterID_W
writeCSV(Cluster_Coord_W,d(:,3),Weights{1},path)

% Density
W=round(Weights{2}.*10);
u = repelem(Coord,W,1);
[ClusterID_W, Cluster_Coord_W ] = kmeans(u,Num_Clusters);
path=fullfile(pwd,'ATL_map_files',strcat('Gwinett_WT_',WType{2},'_Centroids_N',num2str(Num_Clusters),'.csv'));
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
    
    Center_X = 751415.9826674;
    Center_Y = 3754662.1551895;
    %Center=[Center_X Center_Y];
    
    X=[Center_X; X];
    Y=[Center_Y; Y];
    Weight=[0; Weight];
    
    T = table(X,Y,Weight);
    writetable(T,path,'WriteRowNames',true) 
end
