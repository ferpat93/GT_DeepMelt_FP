% Statistical Analysis
close all
%function []=StatSM()

% Number of Folders
nFolders=1;
fi=0;

%Parameters
Treatments={'Glucose_200mM_NaCl_200mM','Glucose_100mM_NaCl_200mM','Glucose_100mM','Glucose_200mM'};
name_file='Evolution_over_time';
ResultsFolder=['C:\Users\lfp3\Desktop' filesep 'StatSM-TotalsToulouse'];
mkdir(ResultsFolder); %Create Results Folder

SM=struct;
BLOB=struct;
AGAR=struct;
for t = 1 : length(Treatments)
    eval(strcat('SM.',Treatments{t},'=[];'))
    eval(strcat('BLOB.',Treatments{t},'=[];'))
    eval(strcat('AGAR.',Treatments{t},'=[];'))
end

folder=cell(nFolders,1);

% Folder F
for nF=1:nFolders
    
    % Get a list of all files and folders in this folder.
    folder{nF}=uigetdir(pwd, strcat('Specify folder of results number: ',num2str(nF)));
    
end

for nF=1:nFolders
    
    % Get a list of all files and folders in this folder.
    files = dir(folder{nF});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    % Print folder names to command window.

    c=struct2cell(subFolders);
    
    names=c(1,3:end);
    UnusedFolders=ones(1,length(names));
    
    for t = 1 : length(Treatments)
        replicates=find(UnusedFolders.*contains(names,Treatments{t}));
        
        for ir=1:length(replicates)
            replicate=replicates(ir);
            UnusedFolders(replicate)=0;
            
            fullpath=fullfile(folder{nF},subFolders(replicate+2).name);
            
            if (exist(fullfile(fullpath,strcat(subFolders(replicate+2).name(9:end),'-',name_file,'.mat')),'file'))>0
                route=fullfile(fullpath,strcat(subFolders(replicate+2).name(9:end),'-',name_file));
            else
                route=fullfile(fullpath,name_file);
            end
            
            str = load(route);
            data_r=str.entities;
            
            %STORE SM
            data_i_SM=NaN(1,700);
            data_i_SM(1:length(data_r(1,:)))=data_r(1,:);
            eval(strcat('SM.',Treatments{t},'(end+1,:)=data_i_SM;'))
            
            %STORE BLOB
            data_i_BLOB=NaN(1,700);
            data_i_BLOB(1:length(data_r(1,:)))=data_r(3,:);
            eval(strcat('BLOB.',Treatments{t},'(end+1,:)=data_i_BLOB;'))
            
            %STORE AGAR
            data_i_AGAR=NaN(1,700);
            data_i_AGAR(1:length(data_r(1,:)))=data_r(2,:);
            eval(strcat('AGAR.',Treatments{t},'(end+1,:)=data_i_AGAR;'))
        end       
    end
end


Entity={'SM','BLOB','AGAR'};
figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
for E=1:length(Entity) %Choose SM or BLOB    
    for t = 1:length(Treatments)
        eval(strcat('Data=',Entity{E},'.',Treatments{t},';'))  
        nRep=length(Data(:,1));
        %disp(strcat(Treatments{t},'  number of rep:  ',num2str(nRep)));
        fi=fi+1;
        subplot(1,3,E)
        hold on
        Mean=nanmean(Data(:,1:450),1);
        Std=nanstd(Data(:,1:450),1,2); 
        plot(Mean);
        xlim([0 450]);
        ylim([0 100]);
    end
end

legend(Treatments)
%Save figures to desktop
saveas(gcf,fullfile(ResultsFolder,'StatisticAnalysis-homogeneous.tif'),'tiffn')