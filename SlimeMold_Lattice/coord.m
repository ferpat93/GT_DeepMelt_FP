function [x1,y1,x2,y2]=coord(n,NX,NY,D)


[i1,j1,i2,j2]=source(n,NX,NY); % gets id of starting and ending node

DX=cos(deg2rad(30))*D; % X distance between nodes
DY=sin(deg2rad(30))*D; % Vertical change in coordinate (odd nodes)

%POINT 1

x1=(j1-1)*DX; % Sets X coordinate for point 1

if mod(j1,2)==1 % If n is odd
    y1=(i1-1)*D+DY;
else
    y1=(i1-1)*D;% If n is even
end

%POINT 2

x2=(j2-1)*DX; % Sets X coordinate for point 2

if mod(j2,2)==1 % If n is odd
    y2=(i2-1)*D+DY;
else
    y2=(i2-1)*D;% if n is even
end