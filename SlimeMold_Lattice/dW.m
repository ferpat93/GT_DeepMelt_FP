function [topology]=dW(matrix,NX,NY,D)

topology=zeros(length(matrix(:,1)),3);

for j=1:length(matrix(:,1)) % Runs on every edge
    [x1,y1,x2,y2]=coord(j,NX,NY,D); %Gets initial and final point coordinates
    topology(j,:)=[(x1+x2)*0.5 (y1+y2)*0.5 matrix(j,1)];
end

