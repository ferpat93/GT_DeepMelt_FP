%% Complimentary Spatial Analysis Launch
clear;
clc;

Treatments={'Glucose_200mM_NaCl_200mM','Glucose_100mM_NaCl_200mM','Glucose_100mM','Glucose_200mM'};
DistanceArray=cell(4,1);
% Get a list of all files and folders in this folder.
ResultsFolder='C:\Users\lfp3\OneDrive - Georgia Institute of Technology\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS';

files = dir(ResultsFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.

%LOOP to Classify and Analyze

for k = 3 : length(subFolders)
    
	fprintf('Sub folder #%d = %s\n', k-2, subFolders(k).name);
    Name=subFolders(k).name;
    Specific_ResultsPath=[ResultsFolder filesep subFolders(k).name];
    route=fullfile(Specific_ResultsPath,strcat(Name(9:end),'-Spatial_analysis.mat'));
    if (exist(route,'file'))>0
        disp('exists')
        str = load(route);
        DA=str.DistAngle(1,:);
        newrow=zeros(1,1000);
        newrow(1:length(DA))=DA;
        for t=1:4
            if contains(Name,Treatments{t})
                DistanceArray{t}=[DistanceArray{t};newrow];
                break
            end
        end

    end

  
end

figure (30); % Spatial plots
hold on
% Maximize the figure. 
set(gcf, 'Position', get(0, 'ScreenSize')); 
set(gcf,'name','Distance to food vs time') 

for t=1:4
    subplot(2,2,t)
    title(Treatments{t});
    hold on
    for i=1:length(DistanceArray{t}(:,1))
        lasti=find(DistanceArray{t}(i,:),1,'last');
        plot(DistanceArray{t}(i,1:lasti))    
    end
    xlabel('Frame')
    ylabel('Distance (pixels)')
    ylim([0 450])
    xlim([0 300])
end

saveas(gcf,fullfile(ResultsFolder,strcat(Name,'-DistanceToFood.tif')))
