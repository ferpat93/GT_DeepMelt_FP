% Process Experiment

%[nc,ma]=PaceParalleltoolbox_r2016b();

DataFolder = 'D:\Slime_Mold_Network\Data\Fusion';
ResultsFolder = 'D:\Slime_Mold_Network\Results\Fusion';
ROI_filename = 'ROI_file.csv';
ImgType = '*.jpg';

Cases={'Control','Control','Glucose','Glucose','NaCl','NaCl'};

Folders=GetSubfolders(DataFolder);

for f=3:numel(Folders) % For each Folder of images (each may contain several dishes)
    
    disp(strcat('Started processing folder: ',Folders{f}))
    
    % Check if ROI's have been found
    if (exist(fullfile(DataFolder,Folders{f},ROI_filename))==2) % If exist do
        
        ROI = csvread(fullfile(DataFolder,Folders{f},ROI_filename));      
        nE = size(ROI,1);
        srcFiles = dir(fullfile(DataFolder,Folders{f},ImgType));
        nI=numel(srcFiles);
        [imageIndexes,diff]=sortTime(srcFiles);
        nFrames = numel(imageIndexes);
        
        histogram(diff)
        drawnow()
        
        for e=2:3%nE % Loop over experiments of the folder
            
            disp(strcat('Started processing experiment #',num2str(e)))            
            tic
            
            centroids=[];
            mask = createCircleMask(ROI(e,3));
                             
            Clusters_labeled = cell(nFrames,1);
            VeinsDist = cell(nFrames,1);
            Graphs = cell(nFrames,1);
            
            % parfor or for (images)
            i=1;
            Image = imread(fullfile(DataFolder,Folders{f},srcFiles(imageIndexes(i)).name));
            Image = Image(ROI(e,2):ROI(e,2)+ROI(e,3)-1,ROI(e,1):ROI(e,1)+ROI(e,3)-1,:);
            [Clusters_labeled{i},VeinsDist{i},Graphs{i},centroids] = getNetwork_Fusion(i,Image,mask,centroids);
%                 figure
%                 plotG(Graphs{i},Clusters_labeled{i})
            
            parfor i=2:nFrames % (loop over images - time)
                Image = imread(fullfile(DataFolder,Folders{f},srcFiles(imageIndexes(i)).name));
                Image = Image(ROI(e,2):ROI(e,2)+ROI(e,3)-1,ROI(e,1):ROI(e,1)+ROI(e,3)-1,:);
                [Clusters_labeled{i},VeinsDist{i},Graphs{i},~] = getNetwork_Fusion(i,Image,mask,centroids);
%                 figure
%                 plotG(Graphs{i},Clusters_labeled{i})
            end
            
            save(fullfile(ResultsFolder,strcat(Folders{f},'_E',num2str(e),'_',Cases{e},'.mat')),'Clusters_labeled','VeinsDist','Graphs','diff','-v7.3')
            
            %RGB = Intensity2RGB(VeinsDist,Clusters_labeled,mask);
            %pathVideo = fullfile(ResultsFolder,strcat(Folders{f},'_E',num2str(e),'_',Cases{e},'.avi'));
            %CreateVideoFromCell(RGB,pathVideo,3)
            toc
        end % For Experiments in folder
        
    end 
    
end % For Folders

%PlayVideo(pathVideo)
%}


%% Auxiliary functions

function [names]=GetSubfolders(path)
        d = dir(path);
        isub = [d(:).isdir]; %# returns logical vector
        names = {d(isub).name}';
        names(ismember(names,{'.','..'})) = [];
end
function [] = CreateVideoFromCell(C,path,FrameRate)
    outputVideo = VideoWriter(path);
    outputVideo.FrameRate = FrameRate;
    open(outputVideo)

    for ii = 1:numel(C)
       img = C{ii};
       writeVideo(outputVideo,img)
    end
    
    close(outputVideo)
end
function [] = PlayVideo(path)
    
    Video = VideoReader(path);

    ii = 1;
    while hasFrame(Video)
       mov(ii) = im2frame(readFrame(Video));
       ii = ii+1;
    end

    figure 
    imshow(mov(1).cdata, 'Border', 'tight')
    movie(mov,1,Video.FrameRate)
end
function [mask] = createCircleMask(L)
    
    [xx,yy] = meshgrid(1:L,1:L);
    mask = false(L,L);
    mask = mask | hypot(xx - L/2, yy - L/2) <= 0.99*(L/2);

end
function [RGB] = Intensity2RGB(Veins,Clusters,mask)
    
    nColors =256;
    colors  = jet(nColors);
    colors(1,:) = [1,1,1];
    
    MV = max(cellfun(@(x) max(x(:)),Veins,'UniformOutput',true)); % maximum overall value
    K = (nColors-1)/MV;

    RGB=cell(numel(Veins,1));
    size_I = size(mask);

    
    for i=1:numel(Veins)       
        vi = ind2rgb(round(K.*Veins{i}+1),colors);
        
        vi = reshape(vi,[],3);
        vi(find(Clusters{i}),:) = 125/255;
        vi(find(~mask),:) = 0;     
        vi = reshape(vi,[size_I 3]);       

        RGB{i} = imresize(uint8(round(255.*vi)),0.5);
    end
end
function []=plotG(G,labeled_clusters)
    
    brown=[173,96,10]./255;
    c=copper(4);
    
    figure()
    plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'EdgeColor',brown,'NodeColor',c(2,:),...
        'NodeLabel',{},'LineWidth',G.Edges.Width./3,'MarkerSize',G.Nodes.Size./3);
    %plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'EdgeColor',brown,'NodeColor',c(G.Nodes.Type+1,:),...
    %    'NodeLabel',{},'LineWidth',G.Edges.Width./3,'MarkerSize',G.Nodes.Size./3);
    set(gca, 'YDir','reverse')
    axis equal
    title('Graph')
    L=length(labeled_clusters(:,1));
    xlim([0 L]) 
    ylim([0 L])
    
    hold on
    
    nClusters = max(labeled_clusters(:));
    for c=1:nClusters
        [pgonY, pgonX] = ind2sub(size(labeled_clusters),find(bwperim(labeled_clusters==c)));
        [ip] = boundary(pgonX, pgonY);
        plot(polyshape(pgonX(ip),pgonY(ip),'Simplify',true),'FaceAlpha',0.9)
    end
    
    
end
function [indexes,diff]=sortTime(srcFiles)
    datetomin=6.94444e-4;
    DataFolder=srcFiles(1).folder;
    nFiles=numel(srcFiles);
    dates=cell(nFiles,1);
    for i=1:nFiles
      dates{i}=extractfield(imfinfo(fullfile(DataFolder,srcFiles(i).name)),'DateTime');  
    end
    
    d=cell2mat(cellfun(@(x) datenum(x,'yyyy:mm:dd HH:MM:SS')-2019*365,dates,'UniformOutput',false));
    [~,indexes]=sort(d);
    diff=round((d(indexes(2:end))-d(indexes(1:end-1)))./datetomin);
end
