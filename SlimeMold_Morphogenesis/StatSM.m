% Statistical Analysis

%function []=StatSM()

% Number of Folders
nFolders=3;
fi=0;

%Parameters
Treatments={'Control','Glucose_100','Glucose_200','NaCl_100'};
name_file='Evolution_over_time';
ResultsFolder=['C:\Users\lfp3\Desktop' filesep 'StatSM-figure'];
mkdir(ResultsFolder); %Create Results Folder

SM=struct('Control',[],'Glucose_100',[],'Glucose_200',[],'NaCl_100',[]);
BLOB=struct('Control',[],'Glucose_100',[],'Glucose_200',[],'NaCl_100',[]);
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
    
    names=c(1,4:end);
    
    for t = 1 : length(Treatments)
        replicates=find(contains(names,Treatments{t})); 
        for ir=1:length(replicates)
            replicate=replicates(ir);
            fullpath=fullfile(folder{nF},subFolders(replicate+2).name);
            
            if (exist(fullfile(fullpath,strcat(subFolders(replicate+2).name(9:end),'-',name_file,'.mat')),'file'))>0
                route=fullfile(fullpath,strcat(subFolders(replicate+2).name(9:end),'-',name_file));
            else
                route=fullfile(fullpath,name_file);
            end
            
            str = load(route);
         
            [data_r_SM,data_r_BLOB]=SMFIT(str.entities');
            
            %STORE SM
            data_i_SM=NaN(500,5);
            data_i_SM(1:length(data_r_SM(:,1)),:)=data_r_SM;
            eval(strcat('SM.',Treatments{t},'(:,end+1,:)=data_i_SM;'))
            
            %STORE BLOB
            data_i_BLOB=NaN(500,5);
            data_i_BLOB(1:length(data_r_BLOB(:,1)),:)=data_r_BLOB;
            eval(strcat('BLOB.',Treatments{t},'(:,end+1,:)=data_i_BLOB;'))            
        end       
    end
end


Entity={'SM','BLOB'};
columns={'Data','Polynomial Fit','Fourier Fit','Poly Derivative','Fourier Derivative'};
x=linspace(1,500,500)';
tot=cell(2,2,length(Treatments));

for E=1:length(Entity) %Choose SM or BLOB
    r=(E-1)*length(Entity)+1;
    
    for t = 1:length(Treatments)
        eval(strcat('Data=',Entity{E},'.',Treatments{t},';'))  
        nRep=length(Data(1,:,1));
        
        fi=fi+1;
        figure(fi)
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
        for sp=1:5
            dim = [(.1625*sp) .93 .1 .05];
            annotation('textbox',dim,'String',columns{sp},'FitBoxToText','on')
            
            subplot(2,5,sp)
            hold on
            for r=1:nRep     
                plot(x,Data(:,r,sp))
            end
            
            Mean=nanmean(Data(:,:,sp),2);
            Std=nanstd(Data(:,:,sp),1,2); 
            
            y=reshape(Data(:,:,sp),[],1);
            xx=repmat(x,nRep,1);
            D=[xx y];
            D(any(isnan(D), 2), :) = [];
            
            Ffitresult = fit(D(:,1),D(:,2),'fourier8');
            Pfitresult = fit(D(:,1),D(:,2),'poly9');

            
            subplot(2,5,sp+5)
            hold on
            plot(x,Mean)
            plot(Ffitresult,'b-','predobs')
            plot(Pfitresult,'r-','predobs')
            legend('off')
   
            if or(sp==1,sp==5)
                ift=round(sqrt(sp));
                tot{E,ift,t}=Ffitresult;
            end
                
        end
        
        legend({'Mean','Fourier Series','Polynomial','Fourier ConfBounds','Poly ConfBounds'},'Position',[0 0.25 0.1 0.05])
        set(legend,'FontSize',9)
        dim = [0.05 0.75 0.1 0.05];
        annotation('textbox',dim,'String',{Entity{E},Treatments{t}},'FitBoxToText','on')
        
        saveas(gcf,fullfile(ResultsFolder,strcat(Entity{E},'-',Treatments{t},'-Stats.tif')))
        
    end
end

figure(9)
set(gcf,'units','normalized','outerposition',[0 0 1 1])

for E=1:length(Entity)
    for c=1:2
        for t = 1:length(Treatments)
            subplot(length(Entity),2,(E-1)*2+c)
            hold on
            y=tot{E,c,t}(x);
            plot(x,y);
        end
    end
end

legend(Treatments,'Position',[0 0.5 0.1 0.05])
annotation('textbox',[0.05 0.25 0.1 0.05],'String','Residuum','FitBoxToText','on')
annotation('textbox',[0.05 0.75 0.1 0.05],'String','Slime Mold','FitBoxToText','on')
annotation('textbox',[0.25 0.92 0.1 0.05],'String','Percentage of total area','FitBoxToText','on')
annotation('textbox',[0.65 0.92 0.1 0.05],'String','Change over time','FitBoxToText','on')
subplot(2,2,2)
ylim([-0.4 0.4])

%Save figures to desktop
saveas(gcf,fullfile(ResultsFolder,'StatisticAnalysis-Glucose-Salt.tif'),'tiffn')