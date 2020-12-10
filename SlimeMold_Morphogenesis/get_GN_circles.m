function [Circle]=get_GN_circles(ImagesFolder,imgType,type)

srcFiles =dir([ImagesFolder '/' imgType]); % Gets info of the files inside the given folder with the given extension

% Ask for mask 
filename1 = fullfile(strcat(ImagesFolder,'/',srcFiles(1).name)); %get first image in the folder
Image = imread(filename1); % open the image
Image = Image(:,:,1:3); %Erase transparency dimension
%[rows, columns, ~] = size(Image); %Get size of images

f=figure(40);
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
hI=imshow(Image);
OKMask=uicontrol('Style','pushbutton',...
                  'String','OK','Units', 'normalized',...
                  'Position',[0.9 0.1 0.05 0.03],'Visible','on',...
                  'Callback','uiresume(gcbf)');
drawnow();

message = strcat('Click 3 times around the circle of',type);
text(10, 40, message, 'color', 'r', 'FontSize', 12);
uiwait(msgbox(message));

% Now, have the user draw the mask in the imageby clicking 3 points
click=0;
circle=zeros(3,2);
while click<3
    [xi, yi, but] = ginput(1);      % get a point
    if ~isequal(but, 1)             % stop if not button 1
        break
    end
    click = click + 1;
    circle(click,:) = [xi yi];
end % Once we get here, the user has finished drawing the region.

[Cmask,Rmask]=calcCircle(circle); %gets center and radius of the drawn mask

CircleObject = imellipse(gca, [Cmask(1)-Rmask Cmask(2)-Rmask 2*Rmask 2*Rmask]);

uiwait(gcf); 
Circle = createMask(CircleObject,hI); %Creates mask in the desired location
close(f);  

end


