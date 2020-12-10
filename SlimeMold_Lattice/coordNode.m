function [cn]=coordNode(nodes,NX,NY,D)

cn=zeros(NX*NY,3);

DX=cos(deg2rad(30))*D; % X distance between nodes
DY=sin(deg2rad(30))*D; % Vertical change in coordinate (odd nodes)

k=0;

for i=1:NY
    for j=1:NX
k=k+1; 

cn(k,1)=(j-1)*DX; % Sets X coordinate 
cn(k,2)=(i-1)*D+mod(j,2)*DY;% Sets Y coordinate
cn(k,3)=nodes(i,j,1)+0.0001; %sets R
    end
end