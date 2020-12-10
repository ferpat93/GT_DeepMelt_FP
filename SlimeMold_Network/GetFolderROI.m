% Try find circles

%DataFolder = 'D:\Slime_Mold_Network\Data';
DataFolder = 'D:\Slime_Mold_Network\Data\Fusion';

ROI_filename = 'ROI_file.csv';
ImgType = '*.jpg';
Folders=GetSubfolders(DataFolder);

Overwrite = 0;

for e = 1 : numel(Folders)
    
    if and(isfile(fullfile(DataFolder,Folders{e},ROI_filename)),~Overwrite)
        continue
    else
        GetROI(fullfile(DataFolder,Folders{e}),ROI_filename,ImgType);
    end
    
end


%{
% Load 1st image and get 2 circles.

srcFiles = dir(fullfile(Root,'*.cr2'));  % the folder in which ur images exists
nI=length(srcFiles); % Number of slices

Io=imread(fullfile(Root,srcFiles(1).name));  

f = figure;
imshow(Io)
c = uicontrol(f,'Style','popupmenu');
c.Position = [20 75 60 20];
c.String = {'Celsius','Kelvin','Fahrenheit'};
c.Callback = @selection;


%}

%% Auxiliar Functions

function [names]=GetSubfolders(path)
        d = dir(path);
        isub = [d(:).isdir]; %# returns logical vector
        names = {d(isub).name}';
        names(ismember(names,{'.','..'})) = [];
end

function []=GetROI(path,ROI_filename,ImgType)

    srcFiles = dir(fullfile(path,ImgType));  % the folder in which ur images exists
    Io=imread(fullfile(path,srcFiles(1).name));

    f = figure;
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    imshow(Io);
    
    nI = 1;
     % create the figure with drop-down menu
     fontSize = 14;
    message = sprintf('Select the number of Dishes in the image');
    text(10, 40, message, 'color', 'r', 'FontSize', fontSize);
    popup = uicontrol('Style', 'popup',...
                   'String', {'6','20'},...
                   'Position', [15 90 100 50],...
                   'Value',1,'Callback',@popupCallback);  
    % callback for the drop-down menu
    function popupCallback(obj,event)
         sels = get(popup,'String');
         idx  = get(popup,'Value');
         nI = str2double(sels{idx});
         uiresume(gcf);
         close(f)
    end
 
    uiwait(gcf);
    
    ROI = zeros(nI,3);
    for i=1:nI
        
        [Cmask,Rmask]=getMask(Io);
        ROI(i,:) = [Cmask(1:2)-Rmask 2*Rmask];
    end

    %writematrix(ROI,fullfile(path,ROI_filename))
    csvwrite(fullfile(path,ROI_filename),ROI)
    
end

function [Cmask,Rmask]=getMask(Image)
fontSize = 14;
    
    f=figure(40);
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    hI=imshow(Image);
    OKMask=uicontrol('Style','pushbutton',...
                      'String','OK','Units', 'normalized',...
                      'Position',[0.9 0.1 0.05 0.03],'Visible','on',...
                      'Callback','uiresume(gcbf)');
    drawnow();

    message = sprintf('Control-GLucose-NaCL in that order!!! Click 3 times around the circular contour to apply the mask');
    text(10, 40, message, 'color', 'r', 'FontSize', fontSize);
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
    mask = createMask(CircleObject,hI); %Creates mask in the desired location
    close(f);
       
    Rmask = round(Rmask);
    Cmask = round(Cmask);
end

function [centre, radius] = calcCircle(pt)

pt1 = double(pt(1,:));
pt2 = double(pt(2,:));
pt3 = double(pt(3,:));

epsilon = 0.000000001;

delta_a = pt2 - pt1;
delta_b = pt3 - pt2;

ax_is_0 = abs(delta_a(1)) <= epsilon;
bx_is_0 = abs(delta_b(1)) <= epsilon;

% check whether both lines are vertical - collinear
if ax_is_0 && bx_is_0
    centre = [0 0];
    radius = -1;
    warning([mfilename ':CollinearPoints'],'Points are on a straight line (collinear).');    
    return
end

% make sure delta gradients are not vertical
% swap points to change deltas
if ax_is_0
    tmp = pt2;
    pt2 = pt3;
    pt3 = tmp;
    delta_a = pt2 - pt1;
end
if bx_is_0
    tmp = pt1;
    pt1 = pt2;
    pt2 = tmp;
    delta_b = pt3 - pt2;
end

grad_a = delta_a(2) / delta_a(1);
grad_b = delta_b(2) / delta_b(1);

% check whether the given points are collinear
if abs(grad_a-grad_b) <= epsilon
    centre = [0 0];
    radius = -1;
    warning([mfilename ':CollinearPoints'],'Points are on a straight line (collinear).');    
    return
end

% swap grads and points if grad_a is 0
if abs(grad_a) <= epsilon
    tmp = grad_a;
    grad_a = grad_b;
    grad_b = tmp;
    tmp = pt1;
    pt1 = pt3;
    pt3 = tmp;
end

% calculate centre - where the lines perpendicular to the centre of
% segments a and b intersect.
centre(1) = ( grad_a*grad_b*(pt1(2)-pt3(2)) + grad_b*(pt1(1)+pt2(1)) - grad_a*(pt2(1)+pt3(1)) ) / (2*(grad_b-grad_a));
centre(2) = ((pt1(1)+pt2(1))/2 - centre(1)) / grad_a + (pt1(2)+pt2(2))/2;

% calculate radius
radius = norm(centre - pt1);

end

