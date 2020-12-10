function [I_matrix]=plotSpatial(Glucose,NaCl,Mask)

%Color Code

% 0 - Empty
% 1 - Glucose
% 2 - NaCl
% 3 - Mask

[rows,columns]=size(Glucose);
Colors=[255 255 255; 127 255 212 ; 140 0 26 ; 0 0 0].*(1/255);

if isempty(NaCl)
    Image=Glucose+3.*Mask;  
else
    Image=Glucose+2.*NaCl+3.*Mask;  
end
ImageList=reshape(Image,rows*columns,1);
Color_list(:,1:3)=Colors(ImageList+1,:);
I_matrix=reshape(Color_list,rows,columns,3);

