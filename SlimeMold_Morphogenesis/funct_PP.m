%Post-Processing
function []=funct_PP(Name,ResultsPath)  
warning off all
% Continue with the demo.  Do some initialization stuff.
close all;
fontSize = 14;

struct = load(fullfile(ResultsPath,strcat(Name,'-ClassifiedImages')));

Images_BIN=struct.ImagesArray;
Images_TRIN=Images_BIN;
History=zeros(size(Images_BIN));
GR=zeros(size(Images_BIN));

rows=struct.rows;
columns=struct.columns;
nImages=length(Images_BIN(1,:));

%input = inputdlg('Enter time between taken photos [min]:','Time spacing', [1 50]); %Gets time spacing between photos
%dt = str2double(input{:}); 

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
GR(delta==8)=-1; % Refining
GR(delta==-8)=2; % Secondary growth (superpose)
GR(delta==-3)=1; % Primary Growth

entities=[sum(Images_TRIN==1,1); %SM
          sum(Images_TRIN==2,1); %Agar
          sum(Images_TRIN==3,1)].*(100/npix); % Residumm

GR_count=[sum(GR==-1,1); %Refining
          sum(GR==1,1); %primary growth
          sum(GR==2,1)].*(100/npix); %secondary growth
      
 % Rerun from #29  
%GetVideo_Auto(Images_TRIN,GR,rows,columns,ResultsPath,Name)

%% Update ClassifiedImages
ImagesArray=Images_TRIN;
Cfullpath = fullfile(ResultsPath,strcat(Name,'-ClassifiedImages')); % Location of ClassifiedImages
save(Cfullpath,'ImagesArray','rows','columns','-v7.3');

%%

figure (13); % Proportion plots
hold on
% Maximize the figure. 
set(gcf, 'Position', get(0, 'ScreenSize')); 
set(gcf,'name','SLIME MOLD - Image Analysis') 


subplot(2,2,1)
title('Binary Entity evolution', 'FontSize', fontSize);
bar(transpose([entities(1,:); entities(2,:)+entities(3,:)]),'stacked')
legend('Slime Mold','Agar')
ylim([0 100])
xlim([0 nImages])

subplot(2,2,2)
title('Trinary Entity evolution', 'FontSize', fontSize);
bar(transpose([entities(1,:);entities(3,:);entities(2,:)]),'stacked')
legend('Slime Mold','Residumm','Agar')
ylim([0 100])
xlim([0 nImages])

subplot(2,2,3)
title('Binary Growing or Refining', 'FontSize', fontSize);
hold on
plot(GR_count(2,:)+GR_count(3,:))
plot(GR_count(1,:),'--')
legend('Growing','Refining')
ylim([0 max(GR_count(:))])
xlim([0 nImages])

subplot(2,2,4)
title('Trinary Growing or Refining', 'FontSize', fontSize);
hold on
plot(GR_count(2,:))
plot(GR_count(3,:))
plot(GR_count(1,:),'--')
legend('Primary Growth','Secondary Growth','Refining')
ylim([0 max(GR_count(:))])
xlim([0 nImages])

saveas(gcf,fullfile(ResultsPath,strcat(Name,'-Entities-vs-time.tif')))
%Cluster_Analysis(struct,PathName)
Cfullpath = fullfile(ResultsPath,strcat(Name,'-Evolution_over_time')); % Gets info of the files inside the given folder with the given extension
save(Cfullpath,'entities','GR_count','-v7.3');
