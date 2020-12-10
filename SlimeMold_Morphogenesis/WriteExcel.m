%% WRITE SPOT DATA

StructSpotData=load('C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS\OutputData_SpotAnalysis.mat');
SpotData=StructSpotData.Output;
% xlswrite(FILE,ARRAY,SHEET) writes to the specified worksheet.

Treatments={'Glucose_100mM_NaCl_200mM','Glucose_200mM_NaCl_200mM','Glucose_100mM','Glucose_200mM'};
Features={'SlimeMold','Agar','Residuum','Refining','PrimaryGrowth','SecondaryGrowth','DistToFood','propSM','Eccentricity','Circularity','NumClusters'};
SpotResultsFolder='C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS\';

for T=1:4
    filename=fullfile(SpotResultsFolder,Treatments{T});
    for E=1:3
        ENT=reshape(SpotData{T}(E,:,:),[600,20]);
        xlswrite(filename,ENT,Features{E})  
    end  
end

%% WRITE HOMOGENEOUS DATA

StructHomoData=load('C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS\RESULTS\OutputData_HomoAnalysis.mat');
HomoData=StructHomoData.Output;
% xlswrite(FILE,ARRAY,SHEET) writes to the specified worksheet.

Treatments={'Glucose_100mM','Glucose_200mM','Control','NaCl_100mM'};
Features={'SlimeMold','Agar','Residuum','Refining','PrimaryGrowth','SecondaryGrowth','Eccentricity','Circularity','NumClusters'};
SpotResultsFolder='C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS\RESULTS\';

for T=1:4
    filename=fullfile(SpotResultsFolder,Treatments{T});
    for E=1:3
        ENT=reshape(HomoData{T}(E,:,:),[600,20]);
        xlswrite(filename,ENT,Features{E})  
    end  
end