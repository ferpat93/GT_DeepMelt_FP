% Trainer SM
function []=funct_Trainer_Kmeans(ImagesFolder,ResultsFolder,DefaultColormap,UC)

warning off all

imgType = '*.tif'; % change based on image type
  
% 1. RGB Analysis

% 1.1. Go over folders of analysis to retrieve colors
%[ColorsU,mask]=uniquecolors_wmask(ImagesFolder,imgType,AutoMask);
%UC=unique(ColorsU,'rows');

%Erase Black from list
RGB=UC(sum(UC,2)>10,:); 

% Lab Analyses
lab = rgb2lab(RGB);
ab=lab(:,2:3);

colorMaps=cell(3,1);
colorMaps{1}=RGB;
colorMaps{2}=lab;
colorMaps{3}=ab;

centroids=cell(3,1);

for c=1:length(colorMaps)
    
    [class, cs]=kmeans(double(colorMaps{c}),2,'distance','sqEuclidean','Replicates',5);
    
    if class(end)==1 % Identify which class has the majority of the 
        class=2.*(class.^(-1));
        c2=cs(1,:);
        cs=[cs(2,:);c2];
    end 
    colorMaps{c}=[colorMaps{c} class];
    centroids{c}=cs;
end

% Separate pixels and show them



figure (20) %Figure including mask 
hold on
set(gcf, 'Position', get(0, 'ScreenSize')); % Maximize the figure. 

for colMap=1:length(colorMaps)
    
    SM=RGB(colorMaps{colMap}(:,end)==1,:);
    L=ceil(length(SM(:,1))^(1/2));
    Vsm=[SM; zeros(L^2-length(SM),3)];
    Msm=reshape(Vsm,L,L,[]);
    
    subplot(2,3,colMap)
    imshow(Msm)
    
    AG=RGB(colorMaps{colMap}(:,end)==2,:);
    L=ceil(length(AG(:,1))^(1/2));
    Vag=[AG; zeros(L^2-length(AG),3)];
    Mag=reshape(Vag,L,L,[]);
    
    subplot(2,3,colMap+3)
    imshow(Mag)
    
end                

figure(21)
scatter3(colorMaps{1}(:,1),colorMaps{1}(:,2),colorMaps{1}(:,3),[],colorMaps{1}(:,4));

figure(22)
scatter3(colorMaps{2}(:,1),colorMaps{2}(:,2),colorMaps{2}(:,3),[],colorMaps{2}(:,4));

figure(23)
scatter(colorMaps{3}(:,1),colorMaps{3}(:,2),[],colorMaps{3}(:,3));



% Show examples to user to decide
if DefaultColormap==0
    srcFiles =dir([ImagesFolder '/' imgType]);
    nImages=length(srcFiles); %Number of images inside the folder 
    nDisplay=5; %Number of images to display as example
    Id_Display= randperm(nImages,nDisplay);

    f=figure (10); %Figure including mask 
    hold on
    set(gcf, 'Position', get(0, 'ScreenSize')); % Maximize the figure. 

    for k =1:nDisplay %Loop through every image in the folder
        i=Id_Display(k);
        filename = strcat(ImagesFolder,'/',srcFiles(i).name); %Set the current Image identifier
        I = imread(filename); % Read the image 
        [rows, columns, ~] = size(I);
        Il=reshape(I,rows*columns,3);

        subplot(length(colorMaps)+1,nDisplay,i)
        imshow(I);

        for row=1:length(colorMaps)
            [~, Location] = ismember(Il,colorMaps{1}(:,1:end-1), 'rows'); %Corresponding location of colors
            b=[colorMaps{1}; zeros(1,length(colorMaps{1}(1,:)))];
            %if colorMaps{row}(find(colorMaps{1}(:,1:3)==[255 255 255],1),end)==1
            if colorMaps{row}(end,end)==1
               k=2;
            else
               k=1;
            end

            b(colorMaps{row}(:,end)~=k,1:end-1)=0;
            Location(Location==0)=length(b(:,1));
            I_bin_list=b(Location,1:end-1);
            I_bin_matrix=reshape(I_bin_list,rows,columns,3);

            subplot(length(colorMaps)+1,nDisplay,nDisplay*row+i)
            imshow(I_bin_matrix)

        end                
    end


    OKMask=uicontrol('Style','pushbutton',...
                          'String','OK','Units', 'normalized',...
                          'Position',[0.9 0.1 0.05 0.03],'Visible','on',...
                          'Callback','uiresume(gcbf)');
    uiwait(gcf); 
    %----------- Wait until user hits ok
    close(f); 
    [CM,~] = listdlg('PromptString','Select a file:','SelectionMode','single','ListString',{'RGB';'Lab';'ab'});
else
    CM=DefaultColormap;
end

eval(strcat('[Classifier, Accuracy]=Classifier_',num2str(CM),'(colorMaps{',num2str(CM),'})'));

Cfullpath = fullfile(ResultsFolder,'Classifier');
save(Cfullpath,'Classifier','-v7.3');
 











