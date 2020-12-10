%% Get Values

% PATH TO FILES

    % Cavity  
Cavity_Folder='C:\Users\lfp3\Documents\IS_SensitivityAnalysis\Cases\Cavity\';
    % EP
EP_Folder='C:\Users\lfp3\Documents\IS_SensitivityAnalysis\Cases\EP\';

%Initial Coordinates

InitialCoord_Cav=ImportDataCav(fullfile(Cavity_Folder,'CavityC608.csv'),true);
EP_Coord=ImportDataEp(fullfile(EP_Folder,'EPC1.csv'),true);

FileType = '*.csv'; % change based on image type
srcFilesCav =dir([Cavity_Folder '/' FileType]);
srcFilesEP =dir([EP_Folder '/' FileType]);

nModels=625;

EP_Map=zeros(length(EP_Coord(:,1)),nModels);
Cav_coord=zeros(length(InitialCoord_Cav(:,1)),2,nModels);

for i =1:nModels %Loop through every model inside the folders
    Cav_coord(:,:,i)=ImportDataCav(fullfile(Cavity_Folder,srcFilesCav(i).name));
    EP_Map(:,i)=ImportDataEp(fullfile(EP_Folder,srcFilesEP(i).name));
end

% Cfullpath = fullfile(Cavity_Folder,'Cavity_Array'); % Gets info of the files inside the given folder with the given extension
% save(Cfullpath,'Cav_coord','InitialCoord_Cav','-v7.3');
% 
% Cfullpath = fullfile(EP_Folder,'EP_Array'); % Gets info of the files inside the given folder with the given extension
% save(Cfullpath,'EP_Map','EP_Coord','-v7.3');
