% Function Binary Image to undirected graph

function [G] = BinToGraph(skel_dist,labeled_clusters)
    %rich_skel = skel_dist;
    % Rich skel is a skeleton with distance values (vein thickness)
    skel = skel_dist>0;
    
    
    %% Check For empty graph
    if sum((skel-labeled_clusters)>0,'all')<5  % Just the Cluster! - Create graph with a single node.
        centroid_cluster = extractfield(regionprops(labeled_clusters>0,'centroid'),'Centroid');
        ClusterID = max(labeled_clusters(:));
        EdgeTable = table([1 1],1,1, 1, ClusterID, ...
        'VariableNames',{'EndNodes' 'Width' 'NumPixels' 'Length' 'ClusterID'});
        NodeTable = table(centroid_cluster(1),centroid_cluster(2),1,1,...
            'VariableNames',{'X' 'Y','Type'});
        G = graph(EdgeTable,NodeTable);
        return
%     else
%         perimeter=bwperim(labeled_clusters>0);
%         clusterEdges = find(perimeter+skel>1);       
    end
    
    
    % I is the binary skeleton to morph into a graph
    BP = bwmorph(skel,'branchpoints');
    EP = bwmorph(skel,'endpoints');
    
    I_Nodes = EP + BP;
    
    % Get Dilated Node Image - 3
    SE = strel('square', 3);
    LN = bwlabel(I_Nodes);
    LN3 = imdilate(LN,SE);

    %% NODES : Graph Nodes: [ X Y Type* ]    *:comes later
    CC = bwconncomp(I_Nodes);
    nNodes = max(LN3(:));
    GraphNodes = zeros(nNodes,3);
    for n=1:nNodes
        [a,b] = ind2sub(CC.ImageSize,CC.PixelIdxList{n}(1)); % Weak, grabs first only
        GraphNodes(n,:) = [b, a, 0]; % Assumes branchpoint
    end
    
    %% EDGES: ['EndNodes' 'Width' 'NumPixels' 'Length'* 'ClusterID']
    
    branches = (skel - (LN3>0))>0;
    LB = bwlabel(branches);
    % CC = bwconncomp(branches);
    SED = strel('disk', 2);
    LB2 = imdilate(LB,SED);
    
    nBranches = max(LB2(:));
    
    if nBranches>0       
        Edges=zeros(nBranches,5);
        for b = 1:nBranches

            % Get Parent Nodes
            nb = unique(LN3(LB2==b))';
            nb(nb==0)=[];       
            if numel(nb)>2
                LN3(LB2==b);
                nb = nb(1:2); % Really bad way to do it -> OJO!
            elseif numel(nb)<2
                nb = [nb nb];
            end

            % Set Edge thickness
            Ds = skel_dist(LB==b);
            npix = numel(Ds);
            Ds(Ds==1)=[];
            if isnan(mean(Ds))
                D = 2;
            else
                D = 2*mean(Ds);
            end

            % Set Cluster ID
            CID = max(labeled_clusters(LB==b));
            if CID>0           
                GraphNodes(nb,3)=1;
                D = 1;
            end
            if numel([nb, D,npix,CID])==5
                Edges(b,:) = [nb, D,npix,CID];
            end

        end
    end
    
    Edges(any(Edges(:,1).*Edges(:,2)==0,2),:)=[]; % Erase Null edges
    
    if or(isempty(GraphNodes),isempty(Edges))
        centroid_cluster = extractfield(regionprops(labeled_clusters>0,'centroid'),'Centroid');
        ClusterID = max(labeled_clusters(:));
        EdgeTable = table([1 1],1,1, 1, ClusterID, ...
        'VariableNames',{'EndNodes' 'Width' 'NumPixels' 'Length' 'ClusterID'});
        NodeTable = table(centroid_cluster(1),centroid_cluster(2),1,...
            'VariableNames',{'X' 'Y','Type'});
        G = graph(EdgeTable,NodeTable);
        return
    else
        Xc=GraphNodes(:,1);
        Yc=GraphNodes(:,2);
        lengths=((Xc(Edges(:,1))-Xc(Edges(:,2))).^2+(Yc(Edges(:,1))-Yc(Edges(:,2))).^2).^(1/2);
    end
    
    %% Set Graph  
    EdgeTable = table(Edges(:,1:2),Edges(:,3),Edges(:,4), lengths, Edges(:,5), ...
        'VariableNames',{'EndNodes' 'Width' 'NumPixels' 'Length' 'ClusterID'});
    NodeTable = table(GraphNodes(:,1),GraphNodes(:,2),GraphNodes(:,3),'VariableNames',{'X' 'Y','Type'});
    G = graph(EdgeTable,NodeTable);

    G = rmedge(G, 1:numnodes(G), 1:numnodes(G)); % Remove self loops
    
    % Update Size of nodes

    internal_nodes = find(~G.Nodes.Type);
    type_list=(G.Edges.ClusterID==0);
    w = G.Edges.Width;
    K =1; % Constant for size

    NodeSize = ones(nNodes,1);

    for n=1:numel(internal_nodes)
        Node = internal_nodes(n);
        eid = outedges(G,Node);
        Ws = w(eid).*type_list(eid);
        NodeSize(Node) = K * sum(Ws) ./ sum(Ws~=0);    
    end

    G.Nodes.Size = NodeSize;
    G=rmnode(G,find(degree(G)==0));
    %{
    % Set Clusters Size
    clusterNodes = find(isnan(NodeSize));
    for c =1:nClusters
        NodeSize(clusterNodes(c)) = 1.5*((cluster_stats(c).Area)/3.14)^0.5;
    end

    G.Nodes.Size = NodeSize;
 %}   
    
end


% Auxiliar Functions

function [Edges] = GetEdges(skel,Nodes)

    Edges=[];
    
    for i=1:numel(Nodes(:,1))
        for f=i+1:numel(Nodes(:,1))
            Ii=removeP(i,f,Nodes,skel);
            %figure()
            %imshow(I);
            D1 = bwdistgeodesic(Ii, Nodes(i,2), Nodes(i,1), 'quasi-euclidean');
            D2 = bwdistgeodesic(Ii, Nodes(f,2), Nodes(f,1), 'quasi-euclidean');
            D = D1 + D2;
            D = round(D * 8) / 8;     
            D(isnan(D)) = inf;
            skeleton_path = imregionalmin(D);

            path_length = D(skeleton_path);
            path_length = path_length(1);


            if (~isinf(path_length)) % If there is a finite path connecting -> Create Edge
                Edges=[Edges ; [i f path_length]];
                %figure()
                %imshow(skeleton_path);
            end
        end
    end
end

function [I]=removeP(i1,i2,Nodes,I)
    m=size(I,2);
    Nodes([i1 i2],:)=[]; 
    %Nodes(Nodes(:,3)==3,:)=[];
    %Nodes(:,3)=[];  
    for i=1:length(Nodes(:,1))
        if Nodes(i,3)==2
            I(Nodes(i,1)-1:Nodes(i,1)+1,Nodes(i,2)-1:Nodes(i,2)+1)=0; 
        else
            I(Nodes(i,1),Nodes(i,2))=0; 
        end
         
        
    end
end

