%% WRITE HOMOGENEOUS DATA

StructHomoData=load('C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS\RESULTS\OutputData_HomoAnalysis.mat');
HomoData=StructHomoData.Output;
% xlswrite(FILE,ARRAY,SHEET) writes to the specified worksheet.

Treatments={'Glucose_100mM','Glucose_200mM','Control','NaCl_100mM'};
Features={'SlimeMold','Agar','Residuum','Refining','PrimaryGrowth','SecondaryGrowth','Eccentricity','Solidity','Circularity','NumClusters'};
SpotResultsFolder='C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS\RESULTS\';
filename=fullfile(SpotResultsFolder,'GrapherData');

t=linspace(0,35,420);

% Sheets 1 to 4 - Entities over time

header={'Time [min]','Slime Mold','Agar','Residuum'};

for tr=1:4
    
    meanMat=NaN(420,4);
    meanMat(:,1)=t';

    percMat=NaN(420*2,4);
    percMat(:,1)=[t'; flipud(t')];
     
    ENT=HomoData{tr}(1:3,1:420,1:20);
    
    meanMat(:,2:end)=nanmean(ENT,3)';
    percMat(1:420,2:end)=prctile(ENT,25,3)';
    percMat(421:end,2:end)=flipud(prctile(ENT,75,3)');

    
    sheet = tr;
    xlRange = 'A1';
    A=[header;num2cell(meanMat)];
    xlswrite(filename,A,sheet,xlRange)

    xlRange = 'H1';
    A=[header;num2cell(percMat)];
    xlswrite(filename,A,sheet,xlRange)   
    
end


% Sheets 5 to 8 - Growth types

header={'Time [min]','Refining','PrimaryGrowth','SecondaryGrowth'};

for tr=1:4
    
    meanMat=NaN(420,4);
    meanMat(:,1)=t';

    percMat=NaN(420*2,4);
    percMat(:,1)=[t'; flipud(t')];
     
    ENT=cumsum(HomoData{tr}(4:6,1:420,1:20),2);
    
    meanMat(:,2:end)=nanmean(ENT,3)';
    percMat(1:420,2:end)=prctile(ENT,25,3)';
    percMat(421:end,2:end)=flipud(prctile(ENT,75,3)');

    
    sheet = tr+4;
    xlRange = 'A1';
    A=[header;num2cell(meanMat)];
    xlswrite(filename,A,sheet,xlRange)

    xlRange = 'H1';
    A=[header;num2cell(percMat)];
    xlswrite(filename,A,sheet,xlRange)   
    
end

% Sheets 6 to 10 - Shape Indexes

header=['Time [min]' Treatments];

meanMat=NaN(420,4,4);
percMat=NaN(420*2,4,4);

%percMat(:,1)=[t'; flipud(t')];
%meanMat(:,1)=t';
   
for tr=1:4
    ENT=HomoData{tr}(7:end,1:420,1:20);    
    meanMat(:,:,tr)=nanmean(ENT,3)';
    percMat(1:420,:,tr)=prctile(ENT,25,3)';
    percMat(421:end,:,tr)=flipud(prctile(ENT,75,3)');
end

for si=1:4
    
    sheet = si+8;
    
    xlRange = 'A1';
    A=[header;num2cell([t' permute(meanMat(:,si,:),[1 3 2])])];
    xlswrite(filename,A,sheet,xlRange)

    xlRange = 'H1';
    A=[header;num2cell([[t'; flipud(t')] permute(percMat(:,si,:),[1 3 2])])];
    xlswrite(filename,A,sheet,xlRange)   
    
end

%% SPOT EXPERIMENTS

StructSpotData=load('C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS\OutputData_SpotAnalysis.mat');
SpotData=StructSpotData.Output;
% xlswrite(FILE,ARRAY,SHEET) writes to the specified worksheet.

Treatments={'Glucose_100mM_NaCl_200mM','Glucose_200mM_NaCl_200mM','Glucose_100mM','Glucose_200mM'};
Features={'SlimeMold','Agar','Residuum','Refining','PrimaryGrowth','SecondaryGrowth','DistToFood','propSM','Eccentricity','Solidity','Circularity','NumClusters'};
SpotResultsFolder='C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS\';
filename=fullfile(SpotResultsFolder,'GrapherData');

t=linspace(0,35,420);

% Sheets 1 to 4 - Entities over time

header={'Time [min]','Slime Mold','Agar','Residuum'};

for tr=1:4
    
    meanMat=NaN(420,4);
    meanMat(:,1)=t';

    percMat=NaN(420*2,4);
    percMat(:,1)=[t'; flipud(t')];
     
    ENT=SpotData{tr}(1:3,1:420,1:20);
    
    meanMat(:,2:end)=nanmean(ENT,3)';
    percMat(1:420,2:end)=prctile(ENT,25,3)';
    percMat(421:end,2:end)=flipud(prctile(ENT,75,3)');

    
    sheet = tr;
    xlRange = 'A1';
    A=[header;num2cell(meanMat)];
    xlswrite(filename,A,sheet,xlRange)

    xlRange = 'H1';
    A=[header;num2cell(percMat)];
    xlswrite(filename,A,sheet,xlRange)   
    
end


% Sheets 5 to 8 - Growth types

header={'Time [min]','Refining','PrimaryGrowth','SecondaryGrowth'};

for tr=1:4
    
    meanMat=NaN(420,4);
    meanMat(:,1)=t';

    percMat=NaN(420*2,4);
    percMat(:,1)=[t'; flipud(t')];
     
    ENT=cumsum(SpotData{tr}(4:6,1:420,1:20),2);
    
    meanMat(:,2:end)=nanmean(ENT,3)';
    percMat(1:420,2:end)=prctile(ENT,25,3)';
    percMat(421:end,2:end)=flipud(prctile(ENT,75,3)');

    
    sheet = tr+4;
    xlRange = 'A1';
    A=[header;num2cell(meanMat)];
    xlswrite(filename,A,sheet,xlRange)

    xlRange = 'H1';
    A=[header;num2cell(percMat)];
    xlswrite(filename,A,sheet,xlRange)   
    
end

% Sheets 6 to 10 - Shape Indexes

header=['Time [min]' Treatments];

meanMat=NaN(420,4,4);
percMat=NaN(420*2,4,4);

%percMat(:,1)=[t'; flipud(t')];
%meanMat(:,1)=t';
   
for tr=1:4
    ENT=SpotData{tr}(9:end,1:420,1:20);    
    meanMat(:,:,tr)=nanmean(ENT,3)';
    percMat(1:420,:,tr)=prctile(ENT,25,3)';
    percMat(421:end,:,tr)=flipud(prctile(ENT,75,3)');
end

for si=1:4
    
    sheet = si+8;
    
    xlRange = 'A1';
    A=[header;num2cell([t' permute(meanMat(:,si,:),[1 3 2])])];
    xlswrite(filename,A,sheet,xlRange)

    xlRange = 'H1';
    A=[header;num2cell([[t'; flipud(t')] permute(percMat(:,si,:),[1 3 2])])];
    xlswrite(filename,A,sheet,xlRange)   
    
end
