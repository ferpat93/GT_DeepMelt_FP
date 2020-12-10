% Generate SM Figure
close all
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.065 0.065], [0.2 0.05], [0.075 0.04]); % subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
% if ~make_it_tight,  clear subplot;  end

to_load =0; 

if to_load ~=0
    DataFolder = 'D:\Slime_Mold_Network\Data\Fusion';
    
    Folder = '24_12_19_Chamber_2';
    time = 585; e = 1;
    
    %Folder = '24_12_19_Chamber_1';
    %e=1; time = round(8.483333*60);
    

    ROI_filename = 'ROI_file.csv';
    srcFiles_filename = 'srcFiles.mat';
    ImgType = '*.jpg';


    ROI = csvread(fullfile(DataFolder,Folder,ROI_filename));      
    nE = size(ROI,1);
    

    srcFiles = load(fullfile(DataFolder,Folder,srcFiles_filename));
    srcFiles = extractfield(srcFiles,'srcFiles'); srcFiles = srcFiles{1};

    mask = createCircleMask(ROI(e,3));
    path = fullfile(DataFolder,Folder,srcFiles.name{time});

    centroids = [ 10 500; 300 500];

    [BW,I]= getBinary(path,ROI(e,:),centroids);

    [G,labeled_clusters,filledBW,skeleton]=BF_Area_Indexes(BW);
    save('sampleSM','G','labeled_clusters','filledBW','BW','skeleton','I');
    
else
    load('sampleSM');
end
a=1;

%figure;
subplot(1,3,3)
plotG(G,labeled_clusters)
xlim([0 500]); ylim([150 900])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
title('c) Undirected graph')

mask = uint8(createCircleMask(size(I,1)));
B = mask.*superposedImage(I,BW,filledBW,skeleton,labeled_clusters);
B = B(150:900,1:500,:);
%figure;
subplot(1,3,2)
imshow(B)
title('b) Processed image')

%figure;
subplot(1,3,1)
C = mask.*I;
imshow(C(150:900,1:500,:))
title('a) Raw image')

function [B] = superposedImage(I,SM,Filled,skeleton,clusters_labeled)
    
    border = bwperim(SM);
    clusters = clusters_labeled>0;
    
    colors = [0.9 0.67 0; 0 0.4 0.4; 0.5 0.0 0];
    position = [1 1];

    A = cat(3,3.*skeleton,2.*clusters,SM);
    L = max(A,[],3);
    B = labeloverlay(I,L,'Colormap',colors);
    %B = uint8(insertText(B,position,'Control - An = 4','BoxOpacity',1,'BoxColor','w'));
    
end


function []=plotG(G,labeled_clusters)
    
    brown=[173,96,10]./255;
    c=copper(4);
    
    plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'EdgeColor',brown,'NodeColor',brown,...
        'NodeLabel',{},'LineWidth',G.Edges.Width./2,'MarkerSize',1);
    %plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'EdgeColor',brown,'NodeColor',c(G.Nodes.Type+1,:),...
    %    'NodeLabel',{},'LineWidth',G.Edges.Width./3,'MarkerSize',G.Nodes.Size./3);
    set(gca, 'YDir','reverse')
    axis equal
%     title('Graph')
%     L=length(labeled_clusters(:,1));
%     xlim([0 L]) 
%     ylim([0 L])
    
    hold on
    
    nClusters = max(labeled_clusters(:));
    for c=1:nClusters
        [pgonY, pgonX] = ind2sub(size(labeled_clusters),find(bwperim(labeled_clusters==c)));
        [ip] = boundary(pgonX, pgonY);
        plot(polyshape(pgonX(ip),pgonY(ip),'Simplify',true),'FaceAlpha',0.9)
    end
    
    
end



function [G,labeled_clusters,filledBW,skel_refined] = BF_Area_Indexes(BW) % returns a table


        %% Split BW into Veins/Clusters
        
        % Identify Clusters and Measure Vein Width
        D_dish = 90; %mm
        ppmm = size(BW,1)/D_dish;   
        maxVeinSize = round(1*ppmm)+1; % 4mm meaning 8mm min dimension.
         
        seE = strel('disk', double(maxVeinSize)); % erode and dilate image
        EI = imerode(BW,seE);
        seD = strel('disk', (round(1.5*maxVeinSize)+1));
        DI = imdilate(EI,seD);

        clustersI = bwlabel(DI.*(BW)); % Get Clusters
        clusters_label = unique(clustersI.*EI);
        clustersI = imfill(ismember(clustersI, clusters_label(2:end)),'holes');
        labeled_clusters = bwlabel(clustersI);
        Full_BI = (BW + clustersI)>0; % after filling inside of cluster
        Idist = bwdist(~Full_BI); % Binary distance (without clusters!)
        Veins = Idist.*(~clustersI);

        fullBW = logical((clustersI + Veins>0)>0);
        filledBW = imfill(fullBW,'holes');

        %% Distribution of W and L
        
        skel_refined = bwskel(Full_BI>0,'MinBranchLength',10);
        skel_dist = skel_refined.*Idist;
        
        G = BuildGraph(skel_dist,clustersI,1,[53 513]); %% Call to external function      

end

function [mask] = createCircleMask(L)
    
    [xx,yy] = meshgrid(1:L,1:L);
    mask = false(L,L);
    mask = mask | hypot(xx - L/2, yy - L/2) <= 0.99*(L/2);

end


function [G] = BuildGraph(skel_dist,ClustersI,needGraph,centroid)
    %rich_skel = skel_dist;
    % Rich skel is a skeleton with distance values (vein thickness)
    skel = skel_dist>0;
    px1mm = size(skel,1)/90; % number of pixels to 1 mm
    
    G=[]; % Initialize empty graph
    
    % I is the binary skeleton to morph into a graph
    BP = bwmorph(skel,'branchpoints');
    EP = bwmorph(skel,'endpoints');
    
    I_Nodes = EP + BP;
    nNodes = sum(I_Nodes,'all');
    
    % Get Dilated Node Image - 3
    SE = strel('square', 3);
    LN3 = imdilate(I_Nodes,SE);
 
    %% EDGES: ['Width' 'NumPixels']
    branches = (skel - (LN3>0))>0;
    LB = bwlabel(branches);  
    nEdges = max(LB(:));

    Per = regionprops(LB,'Perimeter');
    endpts = bwmorph(branches,'endpoints').*LB;
      
        Edges=zeros(nEdges,3); % Thickness, Le, Lp
        
        for b = 1:nEdges

            % Set Edge thickness
            Ds = skel_dist(LB==b);
            Ds(Ds==1)=[];
            if isnan(mean(Ds))
                D = 2.25;
            else
                D = 2*mean(Ds);
            end

            [rows,cols] = find(endpts==b);
            
            if numel(rows) ~= 2
              Le = NaN;
              Lp= nan;
            else
              Le = sqrt((rows(1)-rows(2))^2+(cols(1)-cols(2))^2);
              Lp = Per(b).Perimeter/2;
            end
          
            Edges(b,:) = [D Le Lp]./px1mm;

        end
      
        if needGraph==1 
            G=GetGraph(Edges,I_Nodes,LB,ClustersI,centroid);
        end
   
end



%% Auxiliar Functions

function [G]=GetGraph(Edges,I_Nodes,LB,ClustersI,centroid)

    [LN,nNodes] = bwlabel(I_Nodes);
    SE = strel('square', 3);
    LN3 = imdilate(LN,SE);

    SED = strel('disk', 2);
    LB2 = imdilate(LB,SED);
    
    %% NODES : Graph Nodes: [ X Y Type* ]    *:comes later  
    labeled_clusters = int8(bwlabel(ClustersI));
    s = regionprops(labeled_clusters,'centroid');
    [~,ind]=min(pdist2(centroid,cat(1,s.Centroid)));
    ClustersI = ClustersI + (labeled_clusters==ind);

    GraphNodes=zeros(nNodes,3);
    GraphNodes(:,1:2) = table2array(regionprops('table',LN,'centroid'));
    GraphNodes(:,3) = ClustersI(sub2ind(size(LN),round(GraphNodes(:,2)),round(GraphNodes(:,1))));
        
    nEdges = size(Edges,1);
    GraphEdges = [zeros(nEdges,2) Edges];
    
    for b = 1:nEdges
        % Get Parent Nodes
        nb = unique(LN3(LB2==b))';
        nb(nb==0)=[];       
        if numel(nb)>2
            nb = nb(1:2); % Really bad way to do it -> OJO!
        elseif numel(nb)<2
            nb = [nb nb];
        end      
        GraphEdges(b,1:2)=nb;
    end   
    
    %% Set Graph  
    EdgeTable = table(GraphEdges(:,1:2),GraphEdges(:,3),GraphEdges(:,4),GraphEdges(:,5), ...
        'VariableNames',{'EndNodes' 'Width' 'Le' 'Lp'});
    NodeTable = table(GraphNodes(:,1),GraphNodes(:,2),GraphNodes(:,3),'VariableNames',{'X' 'Y','Type'});
    G = graph(EdgeTable,NodeTable);
    G = rmedge(G, 1:numnodes(G), 1:numnodes(G)); % Remove self loops
    G=rmnode(G,find(degree(G)==0));
    
end




function [BW,I] = getBinary(path,ROI,centroids)

%% Load Image
    Io = imread(path);
    Io = Io(ROI(2):ROI(2)+ROI(3)-1,ROI(1):ROI(1)+ROI(3)-1,:);
    mask = createCircleMask(ROI(3));
    
%% Preprocess image -> Raw to Full Binary

    % pixels per mm
    D_dish = 90; %mm
    ppmm = size(Io,1)/D_dish;
    
    %%% Enhance contrast - grayscale - binarize

    % Enhance
    sigma = 0.4;
    alpha = 0.675;
    I = locallapfilt(Io, sigma, alpha); % Laplacian filtering -  

    lab = rgb2lab(I); % to LAB
    N=lab(:,:,1)-lab(:,:,3); % Compound 
    GI=(N-min(N,[],'all'))/(max(N,[],'all')-min(N,[],'all')); % Normalize
    GI = imgaussfilt(GI,2); % Apply gaussian filter to blur
    
    % Binarize
    BI = (~imbinarize(GI,'adaptive','ForegroundPolarity','dark','Sensitivity',0.55).*(mask));
    
    se = strel('disk',25);
    GI2 = imcomplement(imtophat(imcomplement(GI),se));

    wsi=3/100; % Watershed value - values around 3 and 6% -> depends on gaussian filter

    % Skeleton from watershed
    skelImage = (watershed(imhmin(imcomplement(GI2),wsi))==0);

    BW = ((BI+skelImage)>0).*mask; % Full binary Image  
    
%% Split Cells and measure indexes
    
    T = ppmm^2 * 50; %threshold value to discard noise: considers anything larger than 30 mm2
          
    if isempty(centroids) % Get centroids if its the 1st frame
        BL = bwlabel(BW);
        numPixels = regionprops('table',BL,'Area','Centroid');
        idx = find(numPixels.Area>T);  
        NC = numel(idx);       
        if (NC>2)
            warning('More than 2 large regions in first image')
            pauskjkjke
        else
            centroids=numPixels.Centroid(idx,:);
        end
    end 
    
    nCells = size(centroids,1); % 1 or 2 initial cells
    
    %% Split cells in two
    
    % Connected Print %
    %BW = ConnectedPrint(BW,T,nCells);   
    BW = ConnectedPrint2(BW,T,centroids);   
    
    CC = bwconncomp(BW);
    NC=CC.NumObjects;
    
    if NC<=nCells
        
        DC=inf(nCells,NC); %Rows are centroids, Cols are CC's
        BW=zeros(size(BW,1),size(BW,2),NC);
        
        for k=1:NC
            [row,col] = ind2sub(size(I),CC.PixelIdxList{k});
            for cs=1:nCells
                DC(cs,k)=min(min(((row-centroids(cs,2)).^2+(col-centroids(cs,1)).^2).^(1/2)),DC(cs,k));
            end
            c=zeros(size(BI));
            c(CC.PixelIdxList{k}) = 1;
            BW(:,:,k)=c;
        end
    
        [m,In]=min(DC,[],2);
        if max(m)>25
            warning(strcat('Ojo! Far away centroids - image: ',num2str(ind_im)))
        elseif and(NC==2,nCells==2)
            if In(2)-In(1)==0
                warning('Ojo! not one to one relation')   
            end
        end
        
    else
        
        warning(strcat('More CCs than initial cells!',num2str(ind_im))) 
        pause
    end
    
end
function [Connected2]=ConnectedPrint2(BW,T,centroids)
    
    nCells = size(centroids,1);
    centroids = round(centroids);
    
    BL = bwlabel(BW);
    numPixels = regionprops('table',BL,'Area');
    regionsArea = sortrows([ numPixels.Area (1:numel(numPixels))'],'descend');
    NCCs = sum(regionsArea(:,1)>T); % Number of cells
    
    DC=inf(nCells,NCCs); %Rows are centroids, Cols are CC's
        
    for cs=1:nCells 
        CI = zeros(size(BW));
        CI(centroids(cs,2),centroids(cs,1))=1;
        CID = bwdist(CI);
      
        for k=1:NCCs
            A = CID.*(BL==regionsArea(k,2));
            DC(cs,k)=min(A(A > 0),[],'all');
        end
    end
    
    [m,In]=min(DC,[],2);
    
    iCC = unique(In);
    NC = numel(iCC); % number of cells in the current step (<=nCells)
    
    if sum(m>25)>0
        disp(' Far Away from the centroids')
    end  

    Connected = zeros(size(BW));
    for i=1:NC 
        Connected = Connected + (BL==regionsArea(iCC(i),2));
    end
    
    % Add Left over regions
    DistConnected=bwdist(Connected);   
    %RegionsToAdd=regionsArea(find(and(regionsArea(:,1)<T,regionsArea(:,1)>15)),:); % Leftover regions larger than 15px but less than T
    maxDist = 25; % Max dist to add  
   
    %% Add small Areas
    Connected2=Connected;
    regionsArea(iCC,:)=[];
    for i=1:size(regionsArea,1)
        if regionsArea(i,1)>15
            D=DistConnected+bwdist(BL==regionsArea(i,2));  % Distance between Connected and region i
            if min(D,[],'all')/2<=maxDist % If its close enough
                Connected2 = (Connected2 + (BL==regionsArea(i,2)) + (D<=1+min(D,[],'all')))>0;
            end
        end
    end
    
    %% Check the new areas are connected to large ones
    
    [BL,NC2] = bwlabel(Connected2);

    if NC2>NC
        Connected2 = bwareafilt(Connected2,NC);
%         numPixels = regionprops('table',BL,'Area');
%         regionsArea = sortrows([ numPixels.Area (1:numel(numPixels))'],'descend');
%         NC2 = sum(regionsArea(:,1)>T); % Number of cells
%         
%         Connected2 = zeros(size(BW));
%         for i=1:NC2
%             Connected2 = Connected2 + (BL==regionsArea(i,2));
%         end
    end
    
    % (at this point we know NC2<=NC)
    
    if NC>nCells %  trouble -> shouldn't be here
        
        DistConnected=bwdist(Connected2); 
        LL = bwlabel(Connected2);      
        %Connected2 = zeros(size(BW));
        
        for i=nCells+1:NC2
            D=DistConnected + bwdist(LL==i);  % Distance between Connected and region i
            Connected2 = (Connected2 + (LL==i) + (D<=2+min(D,[],'all')))>0;
        end  

    end   
        
end

