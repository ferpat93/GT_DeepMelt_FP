%trial

 
I=Entities(:,:,1);
figure()
imagesc(I)
colormap(flipud(parula(4)))
colorbar


D = bwdist(I==1);
figure()
imagesc(D)
colormap(flipud(hot))
colorbar

G=GR(:,:,1);
figure()
imagesc(G)
colormap(flag(4))
colorbar

BuffMat=(G>0).*D;
BuffV=BuffMat(:);
BuffV(BuffV==0)=[];
BuffV(isoutlier(BuffV,'mean'))=[];
maxD=prctile(BuffV,99.5);

figure()
imagesc(BuffMat,[0 maxD])
colorbar
colormap(flipud(hot))



IB=(I==1);
SE = strel('disk',ceil(double(maxD)));
ID = imdilate(IB,SE);
Buff = (ID-IB).*I ;


figure()
imagesc(IB)
colormap(flag(2))

figure()
imagesc(ID)
colormap(flag(2))

figure()
imagesc(ID-IB)
colormap(flag(2))

figure()
imagesc(Buff)
colormap(parula(4))