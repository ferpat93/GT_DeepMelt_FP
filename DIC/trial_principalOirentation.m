A=PS{13,2};

close all;

divMap = [brewermap(101,'RdBu')]; % Divergent Colormap (- to +)
seqMap = [brewermap(101,'YlGn')]; % Sequential Colormap (0 to inf)

ix = 100/size(A,1);
xc = ix/2:ix:100-ix/2;
yc = xc;

iy = 70/40;%size(A,3);
zc = 1.*(iy/2:iy:70-ix/2);

[Grids{1,1}, Grids{1,2}] = meshgrid(xc,zc); % Front
[Grids{2,1}, Grids{2,2}] = meshgrid(yc,zc); % Side
[Grids{3,1}, Grids{3,2}] = meshgrid(xc,yc); % Top


iv = 1; % Major (1), Intermediate (2) or Minor (3) principal strain
P = cat(4,A(:,:,:,iv*3:iv*3+2).*A(:,:,:,11+iv),abs(A(:,:,:,12)-A(:,:,:,14))); % Principal Strain Array

Slices{1} = flipud(permute(P(:,round(size(P,2)/2),1:40,:),[3 1 4 2])); % Front
Slices{2} = flipud(permute(P(round(size(P,2)/2),:,1:40,:),[3 2 4 1])); % Side  
Slices{3} = P(:,:,round(size(P,3)/2),:); % Top
 
components = [ 1 2 ; 3 2 ; 3 1];
sizes = [100 70 100];

figure;
for i=1:3
    u = Slices{i}(:,:,components(i,1));
    v = Slices{i}(:,:,components(i,2));
    shear = imresize(Slices{i}(:,:,4),[sizes(components(i,2)) sizes(components(i,1))]);
    
    f = subplot(1,3,i);   
    imagesc(shear);
    colormap(f,seqMap);
    axis image ; hold on;
    q=quiver(Grids{i,1},Grids{i,2},u,v,4) ;
    q.ShowArrowHead = 'off'; q.Color = 'black';
    set(gca, 'YDir','reverse')
end
