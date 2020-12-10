% Spot Analysis Data Extraction

% Check Progress of analyses
clear 
close all
clc

%Output={'-Evolution_over_time.mat','-Spatial_analysis.mat','-BiDish.mat','-ShapeAnalysis.mat'};
Treatments={'Glucose_100mM_NaCl_200mM','Glucose_200mM_NaCl_200mM','Glucose_100mM','Glucose_200mM'};

ResultsFolder='C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS\';

files = dir(ResultsFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);

% Cell with different treatments
Output=cell(4,1); % Cell array with different treatments
for f=1:4 % Creates a NaN array on every output.
    Output{f}=NaN(12,600,20);
    CountTreat=[0 0 0 0];
end

jhk

%LOOP to go trough folders
for k = 3 : length(subFolders)
   
    fprintf('Sub folder #%d = %s\n', k-2, subFolders(k).name);
    Name=subFolders(k).name(9:end);
    Specific_ResultsPath=[ResultsFolder filesep subFolders(k).name];

    for t=1:4 % For loop to find to which treatment the folder belongs to
        if contains(Name,Treatments{t})
            CountTreat(t)=CountTreat(t)+1;
            break
        end
    end
    
    % Evolution over time
    ET = load(fullfile(Specific_ResultsPath,strcat(Name,'-evolution_over_time'))); % positions 1-3 for entities, 4-6 for GR_count
    Entities=ET.entities;
    GR=ET.GR_count;
    ET_add=[Entities;GR];
    Output{t}(1:6,1:length(GR(1,:)),CountTreat(t))=ET_add;
    
    % Spatial Analysis
    SA = load(fullfile(Specific_ResultsPath,strcat(Name,'-Spatial_analysis'))); % Row 7
    DistAngle=SA.DistAngle;
    Output{t}(7,1:length(DistAngle(1,:)),CountTreat(t))=DistAngle(1,:);
    
    % BiDish
    BD = load(fullfile(Specific_ResultsPath,strcat(Name,'-BiDish'))); % Row 8
    propSM=BD.propSM';
    Output{t}(8,1:length(propSM(1,:)),CountTreat(t))=propSM(1,:);
    
    % Shape Analysis
    FA = load(fullfile(Specific_ResultsPath,strcat(Name,'-ShapeAnalysis'))); % Rows 9-12
    ShapeAnalysis=FA.ShapeAnalysis;
    Output{t}(9:end,1:length(ShapeAnalysis(1,:)),CountTreat(t))=ShapeAnalysis;
    %Output{t}(end,1:length(ShapeAnalysis(1,:)),CountTreat(t))=ShapeAnalysis(end,:).^(-1);
end

Cfullpath = fullfile(ResultsFolder,'OutputData_SpotAnalysis'); % Gets info of the files inside the given folder with the given extension
save(Cfullpath,'Output','-v7.3');
