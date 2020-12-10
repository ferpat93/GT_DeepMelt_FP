function [aPs]=createAP(NX,NY,edges)

aP=zeros(NX*NY); % Matrix of coefficients alpha

for pij=1:NX*NY  %LOOP THAT CREATES aP MATRIX - alpha values are updated after each flow
        
    i=ceil(pij/NX); %Finds row index (Y direction)
    j=pij-(i-1)*NX; %Finds column index (X direction)     
    edgesnearby=neighbours(i,j,NX,NY);% finds the edges adjacent to the given node
    suma=sum(edges(edgesnearby,2));% Finds the sum of alpha of the adjacent nodes (second column of edges) 
    aP(pij,pij)=suma+1e-12; %Assign the suma value to the diagonal of the matrix
    
    for ng=1:length(edgesnearby) %Goes from the first neighbor of the node to the last
        [i1,j1,i2,j2]=source(edgesnearby(ng),NX,NY); %Gets the id of the neighbor edge and then the indexes of the underlying  nodes

        if i1==i && j1==j %if the first node is pij (the current)  then the neighbor is the second
            IDN=(i2-1)*NX+j2; %ID of the node based on row and column
            aP(pij,IDN)=-edges(edgesnearby(ng),2);  %assigns to the position(pij,Pneighbor) the value alpha of the edge connecting them.
        else %the neighbor node is the first one
            IDN=(i1-1)*NX+j1; %ID of the node based on row and column
            aP(pij,IDN)=-edges(edgesnearby(ng),2);  %assigns to the position(pij,Pneighbor) the value alpha of the edge connecting them.
        end 
        
    end
end

aPs=1.*aP;
%aPs=sparse(aP);
%isPositiveDefinite(aPs)
%display(eig(aPs))