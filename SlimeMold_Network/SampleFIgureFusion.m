
LegendTreatments = {'Neutral','Nutritive','Adverse'};

DataFolder = 'D:\Slime_Mold_Network\Data\Fusion';
%DataFolder = '/Users/lfp3/Dropbox (GaTech)/GT/Spring-20/SMN/Data';

ROI_filename = 'ROI_file.csv';
srcFiles_filename = 'srcFiles.mat';
ImgType = '*.jpg';
    
LedgerName = 'Ledger.xlsx';
T = readtable(fullfile(DataFolder,LedgerName));

Folders = unique(T.Folder);

folders_to_plot = 2;
dishes_to_plot = [1 3 5];

sigma = 0.4;
alpha = 0.675;
    
TW = 180; % Time window (before/after fusion)

for f=1:numel(folders_to_plot) % For each Folder of images (each may contain several dishes)
    
    folder = Folders{folders_to_plot(f)};
    
    TF=T(ismember(T.Folder,folder),:); %Smaller table with just experiments of folder 'f'
    
    ROIs = csvread(fullfile(DataFolder,folder,ROI_filename));      
    nE = size(ROI,1);
    
    srcFiles = load(fullfile(DataFolder,folder,srcFiles_filename));
    srcFiles = extractfield(srcFiles,'srcFiles'); srcFiles = srcFiles{1};
    
    t = srcFiles.t; % times of images

    figure;   
    for ei=1:numel(dishes_to_plot) % Loop over experiments of the folder
        e = dishes_to_plot(ei);
        ft = str2double(TF(ismember(TF.N_Petri_Dish,e),:).Time_fusion{:}); % Fusion Time from table
        
        [~,idx]=min(abs(t-(ft+TW)),[],1); % Get indexes of images at intervals requested (tf:tf+2h)              
        Fi = idx(1); % Fusion TIME + TW
        
        %ROI = ROIs(i*2-1,:); 
        ROI = ROIs(e,:); 
        path = fullfile(DataFolder,folder,srcFiles.name{idx(1)});
        Io = locallapfilt(imread(path), sigma, alpha); % Laplacian filtering -  
        Io = Io(ROI(2):ROI(2)+ROI(3)-1,ROI(1):ROI(1)+ROI(3)-1,:);
        mask = createCircleMask(ROI(3));
        
        mask3 = repelem(mask,1,1,3);
        
        %Im = uint8(mask).*Io;
        Im=Io; Im(~mask3)= 255;
        subplot(1,3,ei)
        imshow(Im);
        title(LegendTreatments{ei});
    end
    
end

% folder = 'D:\Slime_Mold_Network\Data\Fusion\24_12_19_Chamber_3';
% images={'IMG_2859.JPG','IMG_2950.JPG','IMG_3434.JPG'}; 
  
function [mask] = createCircleMask(L)
    
    [xx,yy] = meshgrid(1:L,1:L);
    mask = false(L,L);
    mask = mask | hypot(xx - L/2, yy - L/2) <= 0.99*(L/2);

end