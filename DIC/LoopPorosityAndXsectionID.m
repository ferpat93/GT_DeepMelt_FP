%% % Run Analysis 

%% List of Paths to retrieve Data
%Folders={'1ST-WEEK-','2ND-WEEK-','LAST-DAY-'};
Folders={'2ND-WEEK-'};
root='G:\';
%Thresholds = readmatrix('Porosity_thresholds.xlsx','Range','B1:D13');
Thresholds = xlsread('Porosity_thresholds.xlsx',1,'B1:D13');

hws = 8; % half the value used in DIC.
ns = 1.5*hws;

%% Nested Loops to analyze each scan

for f=1:numel(Folders) %Weeks of analysis
    TestPath=fullfile(root,strcat(Folders{f},'SC'),filesep);
    nameTests=GetSubfolders(TestPath);
    
    for t=4:4%numel(nameTests) % Iterates Through Tests

        disp(nameTests{t});      
        TT=Thresholds(t,:);
       
        StepPath=fullfile(TestPath,nameTests{t},filesep);
        nameSteps=GetSubfolders(StepPath);
        
        for s=1:numel(nameSteps) % Iterates Though Steps of a Test
            disp(nameSteps{s}); 
            ImagesFolder=fullfile(StepPath,nameSteps{s},'SlicesY');
            ResultsFolder=fullfile(root,strcat(Folders{f},'OP2'),nameTests{t},nameSteps{s});
            
            Device=load(fullfile(ResultsFolder,'Device.mat'));
            Soil=load(fullfile(ResultsFolder,'Soil.mat'));
            
            [Stack,XCidx,Centroid]=GetVolume_OtherDirection(ImagesFolder,Soil,TT(3));
            Stack(Soil.SDM)=0;
            [P,vol] = GetPorosityMap(Stack,TT(1:2),ns,hws);
            
        end
        
    end
end

%% Auxiliary Functions

function []=SaveMultiTiff(path,array)
    imwrite(array(:,:,1), path)
    for i=2:length(array(1,1,:))
        imwrite(array(:,:,i), path, 'writemode', 'append')
    end                       
end

function []=SaveMultiTiff2(path,array)
    if size(array,3)>1028
        for im=1:2
            imwrite(array(:,:,1+floor((size(array,3)*0.5*(im-1)))), [path '_' num2str(im) '.tif'],'Compression','none')
            for i=2:size(array,3)/2
                imwrite(array(:,:,i+floor((size(array,3)*0.5*(im-1)))), [path '_' num2str(im) '.tif'],'Compression','none', 'writemode', 'append')
            end 
        end 
    else
        imwrite(array(:,:,1), [path '.tif'],'Compression','none')
        for i=2:length(array(1,1,:))
            imwrite(array(:,:,i), [path '.tif'],'Compression','none', 'writemode', 'append')
        end 
    end
end

function [names]=GetSubfolders(path)
        d = dir(path);
        isub = [d(:).isdir]; %# returns logical vector
        names = {d(isub).name}';
        names(ismember(names,{'.','..'})) = [];
end
