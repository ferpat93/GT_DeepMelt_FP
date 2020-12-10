
test = 1;
close all
view = 1; % front
divMap=[[0 0 0]; brewermap(101,'RdBu')];
seqMap=[brewermap(101,'YlGn');[0 0 0]];
means = zeros(13,1);
for test=1:13
    figure(test)
    ang1 = medfilt2(angles{test,view,1});
    ang2 = medfilt2(angles{test,view,2});
    angdif = ang2-ang1;

    subplot(1,3,1); imagesc(ang1,[0 90]); %colormap(seqMap) 
    subplot(1,3,2); imagesc(ang2,[0 90]); %colormap(seqMap)
    subplot(1,3,3); imagesc(angdif,[-90 90]./2); colormap(divMap)
    mean(test) = nanmean(angdif(:));
end

