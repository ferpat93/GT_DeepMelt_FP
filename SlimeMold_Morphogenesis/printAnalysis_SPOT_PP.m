% Plot print Analysis

load('E:\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS\BufferChoiceData.mat')
Treatments={'Glucose_100mM_NaCl_200mM','Glucose_200mM_NaCl_200mM','Glucose_100mM','Glucose_200mM'};

WF=1;

if WF==1
    excelPath='E:\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS\BufferChoiceData.xls';
    nRows=numel(masterArray(:,:,:,1));
    nCol=6;
    M=zeros(nRows,nCol);
    nR=20;
    nT=size(masterArray,3);

    M(:,1)=reshape(repmat(1:4,[nR*nT 1]),[],1);
    M(:,2)=repmat(reshape(repmat(1:nR,[nT 1]),[],1),4,1);
    M(:,3)=repmat((1:nT)',nR*4,1);

    row=1;
    for treat=1:4
        for rep=1:nR
            for t=1:nT
                M(row,4:9)=masterArray(treat,rep,t,:);
                row=row+1;
            end
        end
    end

    toWrite=[M(:,1:2) M(:,3).*(5/60) M(:,7) M(:,4)./(M(:,4)+M(:,5)) M(:,8)./(M(:,8)+M(:,9))];
    Names={'Treatment','Replicate','Time (h)','GrowthExtent (pix)','ChoicePG','RandomPG'};      
    xlswrite(excelPath,toWrite) 
    
end

%% Plots
ZM=size(masterArray);
%[Treat Replicates Time Values]
Data=zeros(ZM(3),9,ZM(1));
t=(1:419)*5/60;
header={'Time [h]','PG','Agar','MigrationRate','G/Ag'};
grapher_path='E:\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS\GrapherChoiceData.xls';

colors=jet(4);

meanMat=NaN(419,5,4);
percMat=NaN(419*2,5,4);

for treat=1:4
    
    M=permute(masterArray(treat,:,:,:),[3 4 2 1]);
    % Arrangement time X data X replicates
    
    %CHOICE
    V=[100.*M(:,1,:)./(M(:,1,:)+M(:,2,:)) 100.*M(:,5,:)./(M(:,5,:)+M(:,6,:)) 20.*M(:,4,:) 100.*(M(:,1,:)+M(:,2,:))./(M(:,5,:)+M(:,6,:))];
    V(1,:,:)=[];
    
    meanMat(:,1,treat)=t';   
    percMat(:,1,treat)=[t'; flipud(t')];
     
    meanMat(:,2:end,treat)=[smooth(nanmean(V(:,1,:),3)) smooth(nanmean(V(:,2,:),3)) smooth(nanmean(V(:,3,:),3)) smooth(nanmean(V(:,4,:),3))];
    percMat(:,2,treat)=[smooth(prctile(V(:,1,:),25,3));smooth(flipud(prctile(V(:,1,:),75,3)))];
    percMat(:,3,treat)=[smooth(prctile(V(:,2,:),25,3));smooth(flipud(prctile(V(:,2,:),75,3)))];
    percMat(:,4,treat)=[smooth(prctile(V(:,3,:),25,3));smooth(flipud(prctile(V(:,3,:),75,3)))];   
    percMat(:,5,treat)=[smooth(prctile(V(:,4,:),25,3));smooth(flipud(prctile(V(:,4,:),75,3)))];
    
    %{
    meanMat(1,:,:)=[];
    meanMat(end,:,:)=[];
    percMat(1,:,:)=[];
    percMat(end,:,:)=[];
    %}
    
    sheet = treat;
    xlRange = 'A1';
    A=[header;num2cell(meanMat(:,:,treat))];
    xlswrite(grapher_path,A,sheet,xlRange)

    xlRange = 'H1';
    A=[header;num2cell(percMat(:,:,treat))];
    xlswrite(grapher_path,A,sheet,xlRange) 
    
end



%% PLOTS 

for treat=1:4
    
    figure(46)
    subplot(2,2,treat)
    hold on
    plot(meanMat(2:end,1,treat),meanMat(2:end,2,treat),'-b')
    plot(meanMat(2:end,1,treat),meanMat(2:end,3,treat),'-r')
    legend('Choice','Random')
    title(Treatments{treat})
    
    plot(percMat(2:end,1,treat),percMat(2:end,2,treat),':b')
    plot(percMat(2:end,1,treat),percMat(2:end,3,treat),':r')  

    % legend('Choice','1st Quart Choice','3rd Quart Choice','Random','3rd Quart Random','3rd Quart Random')
       
    figure (10)
    hold on
    plot(meanMat(2:end,1,treat),meanMat(2:end,4,treat),'-','color',colors(treat,:))
    plot(percMat(2:end,1,treat),percMat(2:end,4,treat),':','color',colors(treat,:))
    
    figure (11)
    hold on
    plot(meanMat(2:end,1,treat),meanMat(2:end,5,treat),'-','color',colors(treat,:))
    plot(percMat(2:end,1,treat),percMat(2:end,5,treat),':','color',colors(treat,:))
    
end

figure (10)
%legend(Treatments)
title('Extent of growth region')

figure (11)
%legend(Treatments)
title('Solidity')

