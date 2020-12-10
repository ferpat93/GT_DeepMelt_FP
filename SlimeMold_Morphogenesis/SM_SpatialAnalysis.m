function [DistToFood,AngleToFood,ClosestPoint,Image]=SM_SpatialAnalysis(circleGlucose,circleNaCl,ResultsPath,Name)
%Define folder
if isempty(circleNaCl)
    nCirc=1;
else
    nCirc=2;
end

%% 

%Import array and GET Growing-Refining

struct = load(fullfile(ResultsPath,strcat(Name,'-ClassifiedImages')));

Images_BIN=struct.ImagesArray;
Images_TRIN=Images_BIN;
GR=zeros(size(Images_BIN));

rows=struct.rows;
columns=struct.columns;
nImages=length(Images_BIN(1,:));

mask=reshape(Images_BIN(:,1)==4,rows,columns);
% Proportions

% Value Codes
% SM =1
% Agar = 2
% Blob = 3
% Mask = 4

npix=(rows*columns)-sum(sum(Images_BIN(:,1)==4)); %Num of pixels excluding the mask;

for p=1:rows*columns %Loop through every pixel
    
    firstC=find(Images_BIN(p,:)==1,1,'first');
    
    if or(firstC==columns,isempty(firstC)) % If there is never SM or only grows at the last picture
      continue
    else
      cs=find(Images_BIN(p,firstC+1:end)==2)+firstC; %find columns to change (equal to 2 after first SM)
      Images_TRIN(p,cs)=3;
    end
   
end

%Growth-Refining-Superpose Method
aa=Images_TRIN.^2;
delta=aa-[aa(:,1) aa(:,1:end-1)];
GR(delta==8)=3; % Refining
GR(delta==-8)=2; % Secondary growth (superpose)
GR(delta==-3)=1; % Primary Growth


%% 

%Open Image andGet circles
%[columnsInImage, rowsInImage] = meshgrid(1:columns, 1:rows);

%[cX_G, cy_G, r_G]=getCircle('Glucose');

%circleGlucose = (rowsInImage - cy_G).^2 ...
%    + (columnsInImage - cX_G).^2 <= r_G.^2;


%if nCirc==2
%    [cX_S, cy_S, r_S]=getCircle('NaCL');
%    circleGlucose = (rowsInImage - cy_G).^2 ...
%    + (columnsInImage - cX_G).^2 <= r_G.^2;

%end

%Gets the matrix indexes of the perimeter of the glucose circle
[FoodPerim(:,1), FoodPerim(:,2)]=find(bwperim(circleGlucose)==1);

assignin('base','FoodPerim',FoodPerim)
%% Iterate over time

GR_matrix=reshape(GR,rows,columns,nImages);
ClosestPoint=zeros(nImages,2); %coordinates of the closest pixel
DistToFood=1000.*ones(1,nImages);
AngleToFood=zeros(2,nImages); % First column for actual  angle, second for angle to food
%distThreshold=0.5*r_G;
%assignin('base','GR_matrix',GR_matrix)

for t=2:nImages
    GRpoints_t=[];
    [GRpoints_t(:,1),GRpoints_t(:,2)]=find(GR_matrix(:,:,t));
    % Find closest point from Growing points to Glucose Circle
    [min_distance,matching_coordinates] = calculate_min_distance(GRpoints_t,FoodPerim);
    [ ~ , coupleCoord]= calculate_min_distance(matching_coordinates(:,1:2),ClosestPoint(t-1,:));
    CloseCandidate=coupleCoord(1,1:2);
    TargetCandidate=matching_coordinates(find(ismember(matching_coordinates(:,1:2),CloseCandidate,'rows')),3:4);

    if min_distance<DistToFood(t-1) % compare with previous one - Store distance to food and angle
        DistToFood(t)=min_distance;
        ClosestPoint(t,:)=CloseCandidate;
        % Find its angle to Glucose Circle Centroid/closest point
        
        AngleToFood(2,t)=atand((TargetCandidate(1)-ClosestPoint(t,1))/(TargetCandidate(2)-ClosestPoint(t,2)));
        AngleToFood(1,t)=atand((ClosestPoint(t,1)-ClosestPoint(t-1,1))/(ClosestPoint(t,2)-ClosestPoint(t-1,2)));
        
       % if t>5  
       %     xang=[ones(4,1) ClosestPoint(t-3:t,1) ];
       %     b=xang\ClosestPoint(t-3:t,2);
       %     AngleToFood(1,t)=atand(b(2));
       % else
       %     AngleToFood(1,t)=atand((CloseCandidate(1)-ClosestPoint(t-1,1))/(CloseCandidate(2)-ClosestPoint(t-1,2)));   
       % end
        
    else
        DistToFood(t)=DistToFood(t-1);
        ClosestPoint(t,:)=ClosestPoint(t-1,:);
        AngleToFood(:,t)=AngleToFood(:,t-1);
    end

    
    if DistToFood(t)==0
        break
    end
end


DistToFood(:,1:2)=[DistToFood(:,3) DistToFood(:,3)];
AngleToFood(:,1:2)=[AngleToFood(:,3) AngleToFood(:,3)];
ClosestPoint(1:2,:)=[ClosestPoint(3,:) ; ClosestPoint(3,:)];
DistToFood(:,t+1:end)=[];
AngleToFood(:,t+1:end)=[];
ClosestPoint(t+1:end,:)=[];

Image=plotSpatial(circleGlucose,circleNaCl,mask);

