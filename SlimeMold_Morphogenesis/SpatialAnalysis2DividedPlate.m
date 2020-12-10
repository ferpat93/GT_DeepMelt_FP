% Spatial Analysis  AD Meeting Toulouse

%% Creates the plots of percentage of SM on each side of the plate , saves a picture inside each results folder

%% Complimentary Spatial Analysis Launch
clear;
clc;

Treatments={'Glucose_200mM_NaCl_200mM','Glucose_100mM_NaCl_200mM','Glucose_100mM','Glucose_200mM'};

% Get a list of all files and folders in this folder.
%BaseFolder=uigetdir(pwd, 'Specify the folder containing the PHOTOS and RESULTS folders for the analysis');
%ResultsFolder=fullfile(BaseFolder,'RESULTS');
ResultsFolder='C:\Users\lfp3\OneDrive - Georgia Institute of Technology\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS';
load('C:\Users\lfp3\OneDrive - Georgia Institute of Technology\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS\circles.mat')

files = dir(ResultsFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.

%LOOP to Go through folders

for k = 3 : length(subFolders)
    
	fprintf('Sub folder #%d = %s\n', k-2, subFolders(k).name);
    Name=subFolders(k).name;
    Specific_ResultsPath=[ResultsFolder filesep subFolders(k).name];
 
    %% Check if base spatial Analysis exists
    nameCI=fullfile(Specific_ResultsPath,strcat(Name(9:end),'-ClassifiedImages.mat'));
    close all
    if (exist(nameCI,'file'))>0 
        BiDish(circles{k-2,1},Specific_ResultsPath,Name(9:end));
    else
        disp('No ClassifiedImages.mat file')

    end
    
    
end