% SandBox

%masterArray

ZM=size(masterArray);
%[Treat Replicates Time Values]

header={'Time [min]','Slime Mold','Agar','Residuum'};


for treat=1:4
    M=permute(masterArray(treat,:,:,:),[3 4 2 1]);
    % Arrangement time X data X replicates
    
    %CHOICE
    V=[M(:,1,:)./(M(:,1,:)+M(:,2,:)) M(:,5,:)./(M(:,5,:)+M(:,6,:)) M(:,4,:) (M(:,1,:)+M(:,2,:))/(M(:,5,:)+M(:,6,:))];
    
    meanMat=NaN(420,5,treat);
    meanMat(:,1,treat)=t';

    percMat=NaN(420*2,5,treat);
    percMat(:,1,treat)=[t'; flipud(t')];
     
    meanMat(:,2:end,treat)=[nanmean(V(:,1,:),3) nanmean(V(:,2,:),3) nanmean(V(:,3,:),3) nanmean(V(:,4,:),3)];
    percMat(:,2,treat)=[prctile(V(:,2,:),25,3);flipud(prctile(V(:,2,:),25,3))];
    percMat(:,3,treat)=[prctile(V(:,3,:),25,3);flipud(prctile(V(:,3,:),25,3))];   
    percMat(:,4,treat)=[prctile(V(:,4,:),25,3);flipud(prctile(V(:,4,:),25,3))];
    percMat(:,5,treat)=[prctile(V(:,5,:),25,3);flipud(prctile(V(:,5,:),25,3))];
    
    sheet = treat;
    xlRange = 'A1';
    A=[header;num2cell(meanMat(:,:,treat))];
    xlswrite(grapher_path,A,sheet,xlRange)

    xlRange = 'H1';
    A=[header;num2cell(percMat(:,:,treat))];
    xlswrite(grapher_path,A,sheet,xlRange) 
    
end


