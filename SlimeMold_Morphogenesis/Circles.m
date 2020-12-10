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

%LOOP to get masks
circles=cell(length(subFolders)-2,2);
imgType = '*.tif'; % change based on image type

for k = 3 : length(subFolders)
    Name=subFolders(k).name;
    circles{k-2,1}=get_GN_circles(fullfile(PhotosFolder,subFolders(k).name),imgType,'Glucose');
    
    if contains(subFolders(k).name,'NaCl')
        circles{k-2,2}=get_GN_circles(fullfile(PhotosFolder,subFolders(k).name),imgType,'NaCL');      
    end
end

save([ResultsFolder filesep 'circles'],'circles','-v7.3');