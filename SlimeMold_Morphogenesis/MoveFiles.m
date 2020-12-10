% Check if results are inside photos folder
CurrentNames={'ClassifiedImages.mat','Classifier.mat','Entitiesvid.avi','Entities-vs-time.tif','Evolution_over_time.mat'};
NewNames={'-ClassifiedImages.mat','-Classifier.mat','-Video.avi','-Entities-vs-time.tif','-Evolution_over_time.mat'};

BaseFolder=uigetdir('C:\Users\lfp3\OneDrive - Georgia Institute of Technology\Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS', 'Specify the folder containing the PHOTOS and RESULTS folders for the analysis');
PhotosFolder=fullfile(BaseFolder,'PHOTOS');
ResultsFolder=fullfile(BaseFolder,'RESULTS');
imgType = '*.tif'; % change based on image type

files = dir(PhotosFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);

masksPath = fullfile(ResultsFolder,'Masks'); % fullpath of the masks
   
if (exist(masksPath,'file'))>0 
    load(masksPath);
else
    masks=cell(length(subFolders),1);
end

%LOOP to Classify and Analyze
for k = 3 : length(subFolders)
   
    fprintf('Sub folder #%d = %s\n', k-2, subFolders(k).name);
    Name=subFolders(k).name;
    Specific_PhotosPath=[PhotosFolder filesep subFolders(k).name];
    Specific_ResultsPath=[ResultsFolder filesep 'Results-' subFolders(k).name];
    
    if exist(Specific_ResultsPath,'dir')==0  
        mkdir(Specific_ResultsPath)
        for f=1:numel(CurrentNames)
            NewPath=[Specific_ResultsPath filesep strcat(Name,NewNames{f})];
            CurrentPath=fullfile(Specific_PhotosPath,CurrentNames{f});
            if (exist(CurrentPath,'file'))>0 
                movefile(CurrentPath,NewPath)
            elseif isempty(masks{k})
                masks{k}=getMask_auto(fullfile(PhotosFolder,subFolders(k).name),imgType,0);
            end
        end
    end
    
    
 
end

save(masksPath,'masks','-v7.3');


