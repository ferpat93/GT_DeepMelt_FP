function [mask]=getMask(Image,AutoMask)
	fontSize = 14;
    
    if AutoMask==0

        f=figure(40);
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        hI=imshow(Image);
        OKMask=uicontrol('Style','pushbutton',...
                          'String','OK','Units', 'normalized',...
                          'Position',[0.9 0.1 0.05 0.03],'Visible','on',...
                          'Callback','uiresume(gcbf)');
        drawnow();

        message = sprintf('Click 3 times around the circular contour to apply the mask');
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
    
    else
        
        Rad =600;
        %c=600;
        %b=175;
        size(Image)
        size(imbinarize(Image))
        %I=Image(c-b:c+b,c-b:c+b,:);
        %Ib=imcomplement(imbinarize(I));
        figure
        imshow(rgb2gray(Image))
        figure
        I=imbinarize(rgb2gray(Image),0.5);
        imshow(I)
        %r=regionprops(Ib(:,:,3),{'Centroid'});
        centers = imfindcircles(I,[0.7*Rad 1.3*Rad]);
        %Cx=c-b+r.Centroid(1);
        %Cy=c-b+r.Centroid(2);
        Cx=centers(1,1);
        Cy=centers(1,2);

        [yDim, xDim, ~] = size(Image); %Get size of images
        [xx,yy] = meshgrid(1:yDim,1:xDim);
        mask = false(xDim,yDim);
        mask = mask | hypot(xx - Cx, yy - Cy) <= Rad;
    end
    
