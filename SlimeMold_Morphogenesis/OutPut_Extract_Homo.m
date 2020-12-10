% Spot Analysis Data Extraction

% Check Progress of analyses
clear 
close all
clc

%Output={'-Evolution_over_time.mat','-Spatial_analysis.mat','-BiDish.mat','-ShapeAnalysis.mat'};
Treatments={'Glucose_100mM','Glucose_200mM','Control','NaCl_100mM'};
ResultsFolder='C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS\RESULTS\';

files = dir(ResultsFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);

% Cell with different treatments
Output=cell(4,1); % Cell array with different treatments
CountTreat=[0 0 0 0];

for k = 3 : length(subFolders)
    Name=subFolders(k).name(9:end)
    for t=1:4 % For loop to find to which treatment the folder belongs
        if contains(Name,Treatments{t})
            CountTreat(t)=CountTreat(t)+1;
            break
        end
    end
end



for f=1:4 % Creates a NaN array on every output.
    Output{f}=NaN(10,600,CountTreat(f)); 
end

Counter=zeros(4,1);

%LOOP to go trough folders
for k = 3 : length(subFolders)
   
    fprintf('Sub folder #%d = %s\n', k-2, subFolders(k).name);
    Name=subFolders(k).name(9:end);
    Specific_ResultsPath=[ResultsFolder filesep subFolders(k).name];
 
    for t=1:4 % For loop to find to which treatment the folder belongs
        if contains(Name,Treatments{t})
            disp(t)
            Counter(t)=Counter(t)+1;
            I=Counter(t);    
            break
        end
    end
    
    % Evolution over time
    ET = load(fullfile(Specific_ResultsPath,strcat(Name,'-evolution_over_time'))); % positions 1-3 for entities, 4-6 for GR_count
    Entities=ET.entities;
    GR=ET.GR_count;
    ET_add=[Entities;GR];
    Output{t}(1:6,1:length(GR(1,:)),I)=ET_add;

    % Shape Analysis
    FA = load(fullfile(Specific_ResultsPath,strcat(Name,'-ShapeAnalysis'))); % Rows 7-10
    ShapeAnalysis=FA.ShapeAnalysis;
    Output{t}(7:end,1:length(ShapeAnalysis(1,:)),I)=ShapeAnalysis;
    %Output{t}(end,1:length(ShapeAnalysis(1,:)),CountTreat(t))=ShapeAnalysis(end,:).^(-1);

    
end


Cfullpath = fullfile(ResultsFolder,'OutputData_HomoAnalysis'); % Gets info of the files inside the given folder with the given extension
save(Cfullpath,'Output','-v7.3');
