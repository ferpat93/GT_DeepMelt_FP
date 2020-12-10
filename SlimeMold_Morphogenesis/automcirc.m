

paths={'/Users/lfp3/Dropbox/IMG_2757.tif','/Users/lfp3/Dropbox/IMG_4163.tif'};

for i=1:length(paths)
    Image = imread(fullfile(paths{i}),'tif');

%    [rows, columns, ~] = size(I)
%    mask=getMask(I);

end

Rad =550;
% +/- 175 (each side)

% c ~ 600 600
c=600;
b=175;

I=Image(c-b:c+b,c-b:c+b,:);
Ib=imcomplement(imbinarize(I));

r=regionprops(Ib(:,:,3),{'Centroid'})

figure(1)
hold on
imshow(Ib(:,:,3))
viscircles([r.Centroid(1) r.Centroid(2)],5)


Cx=c-b+r.Centroid(1)
Cy=c-b+r.Centroid(2)

[yDim, xDim, ~] = size(Image); %Get size of images
[xx,yy] = meshgrid(1:yDim,1:xDim);
mask = false(xDim,yDim);
mask = mask | hypot(xx - Cx, yy - Cy) <= Rad;


figure(3)
hold on
imshow(mask)

viscircles([Cx Cy],Rad)
viscircles([Cx Cy],5)

figure(2)
hold on
imshow(Image)

viscircles([Cx Cy],Rad)
viscircles([Cx Cy],5)
