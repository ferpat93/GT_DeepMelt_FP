function []=BiDish(circleGlucose,ResultsPath,Name)

    %Import array and GET Growing-Refining

    struct = load(fullfile(ResultsPath,strcat(Name,'-ClassifiedImages')));

    Images_BIN=struct.ImagesArray;
    rows=struct.rows;
    columns=struct.columns;
    nImages=length(Images_BIN(1,:));
    npix=(rows*columns)-sum(sum(Images_BIN(:,1)==4)); %Num of pixels excluding the mask;

    %% 
    % Proportions

    % Value Codes
    % SM =1
    % Agar = 2
    % Blob = 3
    % Mask = 4

    %Gets center of SM and glucose, to make line
    InitialSM=reshape(Images_BIN(:,4)==1,rows,columns);
    CenterGlucoseS = regionprops(circleGlucose, 'Centroid');
    CenterSMS = regionprops(InitialSM, 'Centroid');
    CenterGlucose = CenterGlucoseS.Centroid;
    CenterSM = CenterSMS.Centroid;
    dy=CenterSM(2)-CenterGlucose(2);
    dx=CenterSM(1)-CenterGlucose(1);
    k=(CenterGlucose(2)*CenterSM(1)-CenterSM(2)*CenterGlucose(1))/dx;
    Ism=reshape(Images_BIN==1,rows,columns,nImages);
    %% Iterate over time
    propSM=zeros(nImages,1);
    for t=1:nImages
        [xSM, ySM]=find(Ism(:,:,t));
        %   d=(x−x1)*(y2−y1)−(y−y1)*(x2−x1)
        %Side=(xSM-CenterSM(1)).*dy-(ySM-CenterSM(2)).*dx;
        Side=k+(dy/dx).*ySM-xSM;
        propSM(t)=100*sum(Side<0)/length(xSM);  
    end
    assignin('base','propSM',propSM)
    figure (13); % Proportion plots
    hold on
    % Maximize the figure. 
    set(gcf, 'Position', get(0, 'ScreenSize')); 
    set(gcf,'name','SLIME MOLD - Image Analysis') 

    title('BiArea SM Growth over time');
    bar(transpose([propSM; 100-propSM]),'stacked')
    %legend('Slime Mold','Agar')
    ylim([0 100])
    xlim([0 nImages])

    saveas(gcf,fullfile(ResultsPath,strcat(Name,'-BiDish.tif')))

    %figure(130)
    %hold on
    %scatter(xSM(Side<0),ySM(Side<0),'*')
    %scatter(xSM(Side>0),ySM(Side>0),'*')
    %plot([CenterGlucose(2);CenterSM(2)],[CenterGlucose(1);CenterSM(1)])

Cfullpath = fullfile(ResultsPath,strcat(Name,'-BiDish')); % Gets info of the files inside the given folder with the given extension
save(Cfullpath,'propSM','-v7.3');

end