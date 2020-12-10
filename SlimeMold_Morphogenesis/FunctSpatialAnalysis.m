function []=FunctSpatialAnalysis(circles,PhotosFolder,Name,Specific_ResultsPath,k)

imgType = '*.tif'; % change based on image type

srcFiles =dir([fullfile(PhotosFolder,Name) '/' imgType]); % Gets info of the files inside the given folder with the given extension
filename1 = fullfile(strcat(fullfile(PhotosFolder,Name),'/',srcFiles(1).name)); %get first image in the folder
Image = imread(filename1); % open the image
Image = Image(:,:,1:3); %Erase transparency dimension
[DistToFood,AngleToFood,ClosestPoint,~]=SM_SpatialAnalysis(circles{k-2,1},circles{k-2,2},Specific_ResultsPath,Name);

figure (30); % Spatial plots
hold on
% Maximize the figure. 
set(gcf, 'Position', get(0, 'ScreenSize')); 
set(gcf, 'Position', [1 1 1900 700]);
set(gcf,'name','SLIME MOLD - Image Analysis') 

subplot(1,3,1)
%axis('square')
title('Distance from Closer SM entity to Glucose');
plot(DistToFood)
xlabel('Frame')
ylabel('Distance (pixels)')

subplot(1,3,2)
title('Comparison between angle of trajectory and angle to glucose');
%axis('square')
hold on
plot(AngleToFood(1,:))
plot(AngleToFood(2,:))
legend('Angle of Growth','Angle to Glucose')
xlabel('Frame')
ylabel('Degree (ang)')

subplot(1,3,3)
title('Spatial distribution of Slime mold')
%axis('square') % Sets equal scale for axis
% Flip the image upside down before showing it
image(Image)
hold on
%plot(FoodPerim(:,2),FoodPerim(:,1))
plot(ClosestPoint(:,2),ClosestPoint(:,1))

saveas(gcf,fullfile(Specific_ResultsPath,strcat(Name,'-Spatial_Analysis.tif')))

DistAngle=[DistToFood;AngleToFood];
Cfullpath = fullfile(Specific_ResultsPath,strcat(Name,'-Spatial_analysis')); % Gets info of the files inside the given folder with the given extension
save(Cfullpath,'ClosestPoint','DistAngle','-v7.3');
