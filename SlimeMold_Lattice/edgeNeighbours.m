function [EdgesToConnect,remedges]=edgeNeighbours(n,NX,NY,remedges)

EdgesToConnect=[];

for k=1:length(n)
    
[i1,j1,i2,j2]=source(n(k),NX,NY);

vn1=neighbours(i1,j1,NX,NY);
vn2=neighbours(i2,j2,NX,NY);

vn=[vn1 vn2];
vn=vn(vn~=n(k));
commonEdges=intersect(vn,remedges);
EdgesToConnect=[EdgesToConnect ; commonEdges];

[~,i] = ismember(remedges, commonEdges);
remedges(i~=0)=[];

end
