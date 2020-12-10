%% Complimentary Spatial Analysis Launch
clear;
clc;

% Get a list of all files and folders in this folder.
BaseFolder=uigetdir(pwd, 'Specify the folder containing the PHOTOS and RESULTS folders for the analysis');
PhotosFolder=fullfile(BaseFolder,'PHOTOS');
ResultsFolder=fullfile(BaseFolder,'RESULTS');

files = dir(PhotosFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.

%idcs   = strfind(Basefolder,filesep);
%ResultsFolder=[Basefolder(1:idcs(end)-1) filesep strcat('Results-',Basefolder(idcs(end)+1:end))];

load('C:\Users\lfp3\OneDrive - Georgia Institute of Technology\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS\circles.mat')
imgType = '*.tif'; % change based on image type

%LOOP to Classify and Analyze

for k = 3 : length(subFolders)
    close all
    Name=subFolders(k).name;
    fprintf('Sub folder #%d = %s\n', k-2, Name);
    
    Specific_ResultsPath=[ResultsFolder filesep strcat('Results-',subFolders(k).name)];

    nameImage=fullfile(Specific_ResultsPath,strcat(Name,'-Spatial_Analysis.tif'));
    existence=(exist(nameImage, 'file') == 2);
    if (existence==1)
     fprintf('Already checked'); 
    else
        
        srcFiles =dir([fullfile(PhotosFolder,subFolders(k).name) '/' imgType]); % Gets info of the files inside the given folder with the given extension
        filename1 = fullfile(strcat(fullfile(PhotosFolder,subFolders(k).name),'/',srcFiles(1).name)); %get first image in the folder
        Image = imread(filename1); % open the image
        Image = Image(:,:,1:3); %Erase transparency dimension

        tic
        [DistToFood,AngleToFood,ClosestPoint,~]=SM_SpatialAnalysis(circles{k-2,1},circles{k-2,2},Specific_ResultsPath,Name);


        figure (30); % Spatial plots
        hold on
        % Maximize the figure. 
        set(gcf, 'Position', get(0, 'ScreenSize')); 
        set(gcf, 'Position', [1 1 1900 700]);
        set(gcf,'name','SLIME MOLD - Image Analysis') 

        subplot(1,3,1)
        %axis('square')
        title('Distance from Closer SM entity to Glucose');
        plot(DistToFood)
        xlabel('Frame')
        ylabel('Distance (pixels)')

        subplot(1,3,2)
        title('Comparison between angle of trajectory and angle to glucose');
        %axis('square')
        hold on
        plot(AngleToFood(1,:))
        plot(AngleToFood(2,:))
        legend('Angle of Growth','Angle to Glucose')
        xlabel('Frame')
        ylabel('Degree (ang)')

        subplot(1,3,3)
        title('Spatial distribution of Slime mold')
        %axis('square') % Sets equal scale for axis
        % Flip the image upside down before showing it
        image(Image)
        hold on
        plot(FoodPerim(:,2),FoodPerim(:,1))
        plot(ClosestPoint(:,2),ClosestPoint(:,1))
       % image(flip(Image, 1));
       % rows=length(Image(:,1,1));
       % hold on
       % plot(ClosestPoint(:,1),rows-ClosestPoint(:,2),'b-*','linewidth',1.5);
       % set(gca,'Ydir','reverse')

        saveas(gcf,fullfile(Specific_ResultsPath,strcat(Name,'-Spatial_Analysis.tif')))
        %Cluster_Analysis(struct,PathName)

        DistAngle=[DistToFood;AngleToFood];
        Cfullpath = fullfile(Specific_ResultsPath,strcat(Name,'-Spatial_analysis')); % Gets info of the files inside the given folder with the given extension
        save(Cfullpath,'ClosestPoint','DistAngle','-v7.3');
        % plot figure with history
        toc
    end
end