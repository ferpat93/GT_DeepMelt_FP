% Copy OUT files to a single folder

Type = 'Disp_';     
path = 'C:\Users\lfp3\Documents\Local_AutoAbaqus';
%path = pwd;
nameFolder = [Type 'OutFiles'];
dest_folder = fullfile(path,nameFolder);
mkdir(path,nameFolder);



Folders = GetSubfolders(path,Type);

for f=1:numel(Folders)
    
    file=fullfile(path,Folders{f},strcat(Folders{f},'.out')); % OUT file

    stat=copyfile(file,dest_folder);
    
    if stat<1
        f
    end     
    
end

%% Auxiliar Functions

function [names]=GetSubfolders(path,string)
        d = dir(path);
        isub = [d(:).isdir]; %# returns logical vector
        names = {d(isub).name}';
        %names(or(ismember(names,{'.','..'}),~contains(names,string))) = [];
        names(~contains(names,string)) = [];
end 