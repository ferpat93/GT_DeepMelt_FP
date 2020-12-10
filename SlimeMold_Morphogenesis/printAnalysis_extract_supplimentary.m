% Neighborhood probability

% Get a list of all files and folders in this folder.
Treatments={'Glucose_100mM','Glucose_200mM','Control','NaCl_100mM'};
root='E:\';
ResultsFolder=fullfile(root,'Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS\RESULTS\');
files = dir(ResultsFolder);
% Extract only those that are directories.
subFolders = files([files.isdir]);
% Print folder names to command window.

% Loop over folders
Counter=zeros(4,1); % Treatment counter
time_range=[2 420];
nT=time_range(2)-time_range(1)+1;

masterArray = zeros([numel(Treatments),20,nT,6],'single');

for k = 23 : 23 %length(subFolders)    
    %disp('Folder: ') k 
    
    Name=subFolders(k).name(9:end);
    Specific_ResultsPath=fullfile(ResultsFolder,subFolders(k).name);
    
    for treat=1:4 % For loop to find to which treatment the folder belongs
        if contains(Name,Treatments{treat})
            Counter(treat)=Counter(treat)+1;
            rep=Counter(treat);    
            break
        end
    end
    
    [Entities,GR]=GetEntitiesGR(Specific_ResultsPath,Name,time_range); % Entities are time shifted already

    % GR Entities 
        % -1 -> Refinement
        % 0 -> No change
        % 1 -> Primary Growth
        % 2 -> Secondary Growth

     % Format: 
    % Format [ nPG nSG nR Buff_size nAg nRe]

    for t=250:250%nT
        tr=t+time_range(1)-1;
        masterArray(treat,rep,tr,:)=getPrintData(Entities(:,:,t),GR(:,:,t));
    end

end

Colors=[206 185 28; 255 255 255 ; 183 219 127; 0 0 0 ].*(1/255);
 
%% Auxiliary Functions

function [ImagesArray, GR]=GetEntitiesGR(path, name,time_range)
    
    images_path =fullfile(path,strcat(name,'-ClassifiedImages'));
    struct = load(images_path);   
    ImagesArray=struct.ImagesArray;
    rows=struct.rows;
    columns=struct.columns;
       
    if numel(unique(ImagesArray(:,200)))<4
        for p=1:rows*columns %Loop through every pixel
            firstC=find(ImagesArray(p,:)==1,1,'first');
            if or(firstC==columns,isempty(firstC)) % If there is never SM or only grows at the last picture
              continue
            else
              cs=find(ImagesArray(p,firstC+1:end)==2)+firstC; %find columns to change (equal to 2 after first SM)
              ImagesArray(p,cs)=3;
            end
        end
        save(path,'ImagesArray','rows','columns','-v7.3');
    end
    
    GR_path =fullfile(path,strcat(name,'-GR_Images'));
    
    try 
        GRstruct = load(GR_path);
        GR=GRstruct.GR;
    catch 
        %Growth-Refining-Superpose Method
        GR=zeros(size(ImagesArray));
        aa=ImagesArray.^2;
        delta=aa-[aa(:,1) aa(:,1:end-1)];
        GR(delta==8)=-1; % Refining
        GR(delta==-8)=2; % Secondary growth (superpose)
        GR(delta==-3)=1; % Primary Growth
        
        save(GR_path,'GR','-v7.3');
    end
   
    % Resize Image and GR arrays
    ImagesArray=reshape(ImagesArray(:,time_range(1)-1:time_range(2)-1),[rows columns (time_range(2)-time_range(1)+1)]);
    GR=reshape(GR(:,time_range(1):time_range(2)),[rows columns (time_range(2)-time_range(1)+1)]);
 
end

function [v]=getPrintData(I,G)
% Format [ nPG nSG nR Buff_size nAg nRe]
IB=(I==1);
D = bwdist(IB);

Colors=[206 185 28; 255 255 255 ; 183 219 127; 0 0 0 ].*(1/255);

figure(10) % Entitites image
ax1=subplot(2,2,1);
imagesc(I)
colormap(ax1,Colors)
minA = min(I(:));
maxA = max(I(:));
legends = {'Slime Mold','Unexplored Surface','Mucus','Outside'};
hold on;
for K = minA : maxA; hidden_h(K-minA+1) = surf(uint8([K K;K K]), 'edgecolor', 'none'); end
hold off
uistack(hidden_h, 'bottom');
legend(hidden_h, legends(minA:maxA) )
title('Entities distribution')

figure(12) % Distance to SM total
ax3=subplot(2,2,3);
dist_mask=D.*(I<4);
image(dist_mask)
colormap(ax3,hot)
colorbar
title('quien sabe')
 
figure(10) % Change 
ax2=subplot(2,2,2);
image(G+2)
colormap(ax2,flag(4))
minA = min(G(:)+2);
maxA = max(G(:)+2);
legends = {'Refinement','No Change','Primary growth','Secondary growth'};
hold on;
for K = minA : maxA; hidden_h(K-minA+1) = surf(uint8([K K;K K]), 'edgecolor', 'none'); end
hold off
uistack(hidden_h, 'bottom');
legend(hidden_h, legends )
title('Growth/refinement between consecutive frames')

figure(10) % Ditance from growth instances
ax2=subplot(2,2,3);
imagesc((G>0).*D) %% distnce to slime mold
colormap(ax2,hot)
colorbar
title('Slime Mold growth distance (in pixels)')

BuffV=reshape((G>0).*D,[],1);
BuffV(BuffV==0)=[];
BuffV(isoutlier(BuffV,'mean'))=[];
maxD=prctile(BuffV,99.5);

nPG=sum(G(:)==1);
nSG=sum(G(:)==2);
nR=sum(G(:)<0);

SE = strel('disk',ceil(double(maxD)));
Buff=(imdilate(IB,SE)-IB).*I ;

figure(10) % Entitites in buffer
ax4=subplot(2,2,4);
Bufi=Buff-1;
Bufi(Bufi<=0)=3;
image(Bufi)
colormap(ax4,Colors(2:end,:))
minA = min(Bufi(:));
maxA = max(Bufi(:));
legends = {'Unexplored Surface','Mucus','Outside'};
hold on;
for K = minA : maxA; hidden_m(K-minA+1) = surf(uint8([K K;K K]), 'edgecolor', 'none'); end
hold off
uistack(hidden_m, 'bottom');
legend(hidden_m, legends)
title('Entities within growth region')

nAg=sum(Buff(:)==2);
nRe=sum(Buff(:)==3);

GD=maxD/size(I,1); % maximum growth relative to perti dish diameter

v=[nPG nSG nR GD nAg nRe];

end


%{

minA = min(A(:));
maxA = max(A(:));
legends = {'first', 'second', 'third', 'fourth', 'fifth', 'sixth'};
hold on;
for K = minA : maxA; hidden_h(K-minA+1) = surf(uint8([K K;K K]), 'edgecolor', 'none'); end
hold off
uistack(hidden_h, 'bottom');
legend(hidden_h, legends(minA:maxA) )

%}