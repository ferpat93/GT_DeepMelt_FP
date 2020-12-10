function [ ClusterEvolution, Frame ] = CluEvo( ActualCC, ClusterEvolution,k,D,nImages)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
New=[vertcat(ActualCC.Info(:).Centroid) vertcat(ActualCC.Info(:).Area)];

if k==1
    ClusterEvolution=NaN(nImages,ActualCC.NumObjects+1,3);
    
    for np=1:ActualCC.NumObjects
        ClusterEvolution(k,np,:)=New(np,:);
    end
    
else

    NewColumn=find(all(isnan(ClusterEvolution(:,:,1))), 1, 'first');
    Existent=find(~isnan(ClusterEvolution(k-1,1:NewColumn-1,3))); %Column Indexes of 
    
    
    [Check, Dist]=rangesearch([ClusterEvolution(k-1,Existent,1)' ClusterEvolution(k-1,Existent,2)'],New(:,1:2),D);

    ExistentID=(horzcat(Check{:}))';
    Distv=(horzcat(Dist{:}))';
    NewID=[];

    for np=1:ActualCC.NumObjects
        l=length(Check{np});
        NewID=[NewID; np.*ones(l,1)];
    end

    ConnMatrix=sortrows([Distv NewID Existent(ExistentID)']);
    
    UsedNew=[];
    UsedExistent=[];

    for i=1:length(ConnMatrix(:,1))
        if ~or(ismember(ConnMatrix(i,2),UsedNew),ismember(ConnMatrix(i,3),UsedExistent)) %If none of them have been used yet
            ClusterEvolution(k,ConnMatrix(i,3),:)=New(ConnMatrix(i,2),:);
            UsedNew=[UsedNew; ConnMatrix(i,2)];
            UsedExistent=[UsedExistent; ConnMatrix(i,3)];
        end
    end

    %Check New Not used - Create New CLuster
    ToCreate=setdiff(linspace(1,length(New(:,1)),length(New(:,1))),UsedNew);
    for i=1:length(ToCreate)
        ClusterEvolution=[ClusterEvolution NaN(nImages,1,3)];
        ClusterEvolution(k,NewColumn,:)=New(ToCreate(i),:);
        NewColumn=NewColumn+1;
    end
    
    
    %Not Necessary to disable unused Clusters (are already NaN)
   
end
   
