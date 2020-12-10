function [i1,j1,i2,j2]=source(n,NX,NY)

if n<= (NX-1)*NY
%%Horizontal fracture
 
i1=ceil(n/(NX-1));
i2=i1;  

j1=n-(i1-1)*(NX-1);
j2=j1+1;

elseif n<=(2*NX*NY-NY-NX)
%% Vertical Fracture        
  
nv=n-NY*(NX-1);
i1=ceil(nv/NX);
i2=i1+1;

j1=nv-NX*(i1-1);
j2=j1;

else
%% Diagonal Fracture

nd=n-(2*NX*NY-NY-NX);
i1=ceil(nd/(NX-1));
i2=i1+1;

no=nd-((NX-1)*(i1-1));

if mod(no,2)==1
    j1=no; 
    j2=j1+1;

else
    j1=no+1; 
    j2=no;
end

end
    