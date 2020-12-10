% Topological analysis

[Name,Path] = uigetfile('*.txt','Chose Topology file'); 
topo=dlmread(strcat(Path,Name));     % Stores as 'nodes' the  matrix

setX=sort(unique(topo(:,1)));
setY=sort(unique(topo(:,2)));

mX=zeros(length(setX),1);
mY=zeros(length(setY),1);

for k=1:length(setX)
    B=topo(:,1)==setX(k);
    mX(k)=sum(topo(B,3))/sum(B);
end

for k=1:length(setY)
    B=topo(:,2)==setY(k);
    mY(k)=sum(topo(B,3))/sum(B);
end
% 

figure(1)
plot(setX,mX);

figure(3)
plot(setY,mY);


