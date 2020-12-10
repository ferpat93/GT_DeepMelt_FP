
%% Get Videos FusionRegion
%% ~ AFTER - FUSION ~ %%

% Constants
TW = 240; % Time window (before/after fusion)
dt = 10; % Time Interval
timesAF = (0:dt:TW)';
A_int = 1:0.5:4; % Area Intervals for part 1 (relative to initial area)
DishTreatments = {'Control','Control','Glucose','Glucose','NaCl','NaCl'};
Treatments = {'Control','Glucose','NaCl'};

DataFolder = 'D:\Slime_Mold_Network\Data\Fusion';
ResultsFolder = 'D:\Slime_Mold_Network\Results\Fusion';
ROI_filename = 'ROI_file.csv';
srcFiles_filename = 'srcFiles.mat';
ImgType = '*.jpg';
Cases={'Control','Control','Glucose','Glucose','NaCl','NaCl'};
LedgerName = 'Ledger.xlsx';
T = readtable(fullfile(DataFolder,LedgerName));

% Gather Data

Folders=GetSubfolders(DataFolder);
exps = [1 3 6];

for f=4:4%numel(Folders)
    
    TF=T(ismember(T.Folder,Folders{f}),:); %Smaller table with just experiments of folder 'f'   
    ROI = csvread(fullfile(DataFolder,Folders{f},ROI_filename));      
    nE = size(ROI,1);
    
    if (isfile(fullfile(DataFolder,Folders{f},srcFiles_filename))) % If exists load it, otherwise compute it.
        srcFiles = load(fullfile(DataFolder,Folders{f},srcFiles_filename));
        srcFiles = extractfield(srcFiles,'srcFiles');
        srcFiles = srcFiles{1};
    else 
        error
    end
    
    t = srcFiles.t; % times of images
      
    if isfile(fullfile(DataFolder,Folders{f},'FusionRegions.mat'))
        load(fullfile(DataFolder,Folders{f},'FusionRegions.mat'));

        for e=1:2:6%:2:6 % just one per treatment
            ft = str2double(TF(ismember(TF.N_Petri_Dish,e),:).Time_fusion{:}); % Fusion Time from table
            [~,idx]=min(abs(t-(ft:dt:ft+TW)),[],1); % Get indexes of images at intervals requested (tf:tf+2h)

            Images = getImages(fullfile(DataFolder,Folders{f}),srcFiles.name(idx),ROI(e,:),FusionRegionsExperiments(e,:),0:dt:TW); 
            GetVideo(Images,fullfile(DataFolder,Folders{f},strcat(Folders{f},'_e',num2str(e))),1)
        end

        
    end

end


function [Images] = getImages(path,names,ROI,FusionRegions,times)
    
    ROA = FusionRegions{3};
    ROIb = bwperim(ROA);
    bounds = FusionRegions{2};
    SM = FusionRegions{1};
    Io = SM(:,:,1);
    Iob = bwperim(Io);
    
    sigma = 0.4;
    alpha = 0.675;

    Images=uint8(zeros([size(Io,1) size(Io,2) 3 numel(names)]));
    colors = [0.2 0.2 0.2; 0.9 0.9 1; 0.9 0.67 0];
    position = [1 1];
    %figure;
    for i=1:numel(names)
        
        I = imread(fullfile(path,names{i}));
        I = I(ROI(2):ROI(2)+ROI(3)-1,ROI(1):ROI(1)+ROI(3)-1,:);
        I = I(bounds(1):bounds(2),bounds(3):bounds(4),:);
        I = locallapfilt(I, sigma, alpha); % Laplacian filtering - 

        A = cat(3,3*SM(:,:,i),2*Io,~ROA);
        L = max(A,[],3);
        B = labeloverlay(I,L,'Colormap',colors);
        B = uint8(insertText(B,position,strcat('Control - t: ',num2str(times(i)),'m'),'BoxOpacity',1,'BoxColor','w'));
        
        %imagesc(B);
        %[X, ~] = frame2im(getframe(gcf));
        
        %if i==1
        %    Images=zeros([size(X,[1 2]) 3 numel(names)]);
        %end
        Images(:,:,:,i) = B;
        
    end 
    
end

function [] = GetVideo(I,path,FrameRate)

    outputVideo = VideoWriter(path);
    outputVideo.FrameRate = FrameRate;
    open(outputVideo)

    for ii = 1:size(I,4)
       img =  I(:,:,:,ii);
       writeVideo(outputVideo,img)
    end
    
    close(outputVideo)
end

function [names]=GetSubfolders(path)
        d = dir(path);
        isub = [d(:).isdir]; %# returns logical vector
        names = {d(isub).name}';
        names(ismember(names,{'.','..'})) = [];
end 