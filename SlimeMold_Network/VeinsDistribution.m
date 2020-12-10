% Function Binary Image to undirected graph

function [W,Wpar,L,Lpar,indexes,G] = VeinsDistribution(skel_dist,ClustersI,Wbins,Lbins,needGraph,centroid)
    %rich_skel = skel_dist;
    % Rich skel is a skeleton with distance values (vein thickness)
    skel = skel_dist>0;
    px1mm = size(skel,1)/90; % number of pixels to 1 mm
    
    G=[]; % Initialize empty graph
    
    %% Check For empty graph
    if sum((int8(skel)-ClustersI)>0,'all')<5  % Just the Cluster! - Create graph with a single node.
        W=zeros(1,numel(Wbins)-1);
        L=zeros(1,numel(Lbins)-1);
        Wpar=zeros(1,2);
        Lpar=zeros(1,2);
        indexes=zeros(1,4);        
        return
    end
       
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
      
    if nEdges>0       
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
        
        indexes = [nNodes nEdges nansum(Edges(:,2))./px1mm nansum(Edges(:,3))./px1mm]; %'nNodes','nEdges','Le','Lp'
        if size(Edges,1)>10
            [W, ~] = histcounts(Edges(:,1),Wbins,'Normalization','probability');
            pd = fitdist(Edges(:,1),'loglogistic');
            Wpar = [pd.mu pd.sigma];
            [L, ~] = histcounts(Edges(:,3),Lbins,'Normalization','probability');
            pd = fitdist(Edges(:,3),'loglogistic');
            Lpar = [pd.mu pd.sigma]; 
        else
            W=zeros(1,numel(Wbins)-1);
            L=zeros(1,numel(Lbins)-1);
            Wpar = [0 0];
            Lpar = [0 0];
        end
        
        if needGraph==1 
            G=GetGraph(Edges,I_Nodes,LB,ClustersI,centroid);
        end
        
    else
        W=zeros(1,numel(Wbins)-1);
        L=zeros(1,numel(Lbins)-1);
        Wpar = [0 0];
        Lpar = [0 0];
        indexes = [0 0 0 0];
        return
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
