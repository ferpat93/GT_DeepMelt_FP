 % find min
 close all
 R=1000;
 
 md = 375;
 angles = 10:0.1:20;

 [Ls,Lt,H]= getLH(angles,R,md);
 
 figure;
 subplot(3,1,1)
 plot(angles,Ls)
 subplot(3,1,2)
 plot(angles,Lt)
 subplot(3,1,3)
 plot(angles,H)
 
 
 function [L]= getL(angles,R,H)
 
 L=zeros(size(angles));
 for i=1:numel(angles)
     ang = deg2rad(angles(i));
     L(i)= R*ang + (1/sin(ang))*(H-R*(1-cos(ang)));
 end
 
 end
 
  function [Ls,Lt,H]= getLH(angles,R,md)
 
 Ls=zeros(size(angles));
 Lt=zeros(size(angles));
 H=zeros(size(angles));

 for i=1:numel(angles)
     ang = deg2rad(angles(i));
     Ls(i) = md/cos(ang) - R*tan(ang);
     Lt(i)= R*ang + Ls(i);
     H(i)= R*(1-cos(ang)) + Ls(i)*sin(ang);
 end
 
  end
 
  function [H]= getH(angles,R,Ls)
 
 H=zeros(size(angles));
 for i=1:numel(angles)
     ang = deg2rad(angles(i));
     H(i)= R*(1-cos(ang)) + Ls*sin(ang);
 end
 
  end
 