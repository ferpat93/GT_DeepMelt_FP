%% Build 3D Volume

function [Volume,XCidx,Centroid]=GetVolume_OtherDirection(slicesPath,SoilStruct,TT,DB)

warningID = 'images:bwfilt:tie'; % Turn off warning of bwlabel - tie at n-th position
warning('off',warningID);

maskCircle = SoilStruct.SCM;
BoxSize = size(maskCircle,1);
Bounds = SoilStruct.SCB;
DB = DB - Bounds(1,1);

srcFiles = dir(strcat(slicesPath,'\*.tif'));  % the folder in which ur images exists
nSlices=length(srcFiles); % Number of slices

%% Load Volume

Volume=zeros([size(maskCircle) nSlices,],'uint16'); % Initialize DevicePrint

for i = 1:nSlices
    Image=imread(strcat(slicesPath,'\',srcFiles(i).name));
    Volume(:,:,i)=Image(Bounds(2,1):Bounds(2,2),Bounds(1,1):Bounds(1,2)).*uint16(maskCircle); % Crop Image and apply mask
end

%% Device Identification

DeviceSlice=round(BoxSize/2); % Number of slice to start looking for device
SED=5; % Square structure element - Identify Device
MinArea=0.8*(3.14/4)*(10/0.070)^2; % min area -> 80% of the circle diameter  
Centroid = [];
XCidx = zeros(round(DB(2)-DB(1)+1),8,'double');

% -> To the front
for i = DeviceSlice : DB(2) 
    [I,Xci,Centroid]=GetCrossSection(Volume(:,i,:),MinArea,SED,TT,Centroid);
    if ~isempty(I)
        Volume(:,i,:)=I;
        XCidx(uint16(i-DB(1)+1),:) = [i Xci];
    else
        disp(i)
    end   
end

% -> To the back
for i = DeviceSlice-1 : -1 : DB(1)  
    [I,Xci,Centroid]=GetCrossSection(Volume(:,i,:),MinArea,SED,TT,Centroid);
    
    if ~isempty(I)
        Volume(:,i,:)=I;
        XCidx(i-DB(1)+1,:) = [i Xci];
    else
        disp(i)
    end   
end

%XCidx(end,:) = [0 DB 0 0 0 0 0];
XCidx = sortrows(XCidx);

end

function [I,indexes,centroid]=GetCrossSection(ImageO,minArea,SED,TT,centroid)

Image = permute(ImageO,[3 1 2]);
IB=(Image<TT); % Apply threshold
seG = strel('disk',3); 
J = ~imdilate(~IB,seG);

LabeledImage = bwlabel(imfill(bwareafilt(J,[0.8*minArea 30*minArea]),'holes'));
s = table2array(regionprops('table',LabeledImage,'centroid'));

if isempty(centroid)    
    [~,is] = min(abs(s(:,1)-(size(Image,2)/2)));
    centroid = round(s(is,:)); % Weak Check
end

[~,is] = min((s(:,1)-centroid(1)).^2+(s(:,2)-centroid(2)).^2);
IsolatedImage = LabeledImage==is;

%% Dilate grains

se2 = strel('disk',SED); % Creates structure element to homogenize image
BW3 = bwareafilt(imopen(IsolatedImage,se2),1); % Erodes Image with Structure Element
BW3modified=imfill(BW3,'holes'); % Fill Holes... Duh!

windowSize = 31;
kernel = ones(windowSize) / windowSize ^ 2;
blurryImage = conv2(single(BW3modified), kernel, 'same');
binaryImage = blurryImage > 0.5; % Rethreshold

se2 = strel('disk',1); % Creates structure element to homogenize image
binaryImage = bwareafilt(imopen(binaryImage,se2),1); % Erodes Image with Structure Element
binaryImage=imfill(binaryImage,'holes'); % Fill Holes... Duh!

%% Rest of algorithm, printed variables

DeviceIndexes=regionprops('table',binaryImage,'FilledArea','MinorAxisLength','MajorAxisLength','Centroid','Solidity','Orientation');

if ~isempty(DeviceIndexes)
    if DeviceIndexes.FilledArea>=minArea
        I = permute(Image.*uint16(~binaryImage),[2 3 1]);
        indexes = table2array(DeviceIndexes).*[0.07 0.07 0.07 0.07 1 (0.07^2) 1]; % Centroid(2), major-minor axis, orientation, area, solidity
    else
        I=[];
        indexes = [];
    end
    
else
    I=[];
    indexes = [];
end

end
