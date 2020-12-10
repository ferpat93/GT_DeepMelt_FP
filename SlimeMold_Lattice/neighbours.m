function [vn]=neighbours(i,j,NX,NY)


Nh=(NX-1)*NY;
Nvh=(2*NX*NY-NY-NX);

%% Horizontal fractures
 
h1=NX*i-i+j-NX;
h2=h1+1;

%% Vertical Fracture        
  
v1=Nh+j+(i-2)*NX;
v2=Nh+j+(i-1)*NX;

%% Diagonal Fracture

d1=Nvh+(j-1)+(i-2+mod(j,2))*(NX-1);
d2=d1+1;

%% CHECK SPECIAL CASES

if j==1    
    h1=0; 
    d1=0;
elseif j==NX   
    h2=0; 
    d2=0;
end

        
if i==1    
    v1=0; 
    if mod(j,2)==0
        d1=0;
        d2=0;
    end
elseif i==NY   
    v2=0; 
    if mod(j,2)==1
        d1=0;
        d2=0;
    end
end

%% CLEAR AND RETURN NEIGHBOURS VECTOR

a=[h1 h2 v1 v2 d1 d2];
vn=a(a~=0);
