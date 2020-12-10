%% Cavity Analysis

% Cavity  
Cavity_Folder='C:\Users\lfp3\Documents\IS_SensitivityAnalysis\Cases\Cavity';
struct = load(fullfile(Cavity_Folder,'Cavity_Array'));
Cavity_points=struct.Cav_coord;
Cavity_points = cat(1, Cavity_points, Cavity_points(2,:,:));
Cavity_points(2,:,:)=[];
nModels=length(Cavity_points(1,1,:));

%Cavity_points(end+1,:,:)=zeros(size(Cavity_points(1,:,:)));


Cav_indexes=zeros(nModels,4);

r0=1;
 
%figure(69)
%axis equal
%hold on
stats=cell(nModels,1);

for i=1:nModels % 20:20% 
    
    figure(i)
    axis equal
    %hold on
    % Fit Cavity Ellipse and get axis
    ellipse_struct=fit_ellipse(Cavity_points(:,1,i),Cavity_points(:,2,i)); % Struct containing ellipse parameters
    C_Major_axis=ellipse_struct.long_radius;
    C_Minor_axis=ellipse_struct.short_radius;

    %Find Indexes
    Cav_indexes(i,1)=100*(((C_Major_axis*C_Minor_axis)/(r0^2))-1);
    Cav_indexes(i,2)=C_Minor_axis/C_Major_axis; 
    
    %plot(Cavity_points(:,1,i),Cavity_points(:,2,i));
    fullcav=([[Cavity_points(:,1,i);flip(Cavity_points(:,1,i));-1.*Cavity_points(:,1,i);flip(-1.*Cavity_points(:,1,i))],[Cavity_points(:,2,i);flip(-1.*Cavity_points(:,2,i));-1.*Cavity_points(:,2,i);flip(Cavity_points(:,2,i))]
                    ]) ;                      
   
   figure(1)
    hold on
    xlim([0 1.2])
    ylim([0 1.2])
    th = 0:pi/100:pi/2;
    r=1;
    xunit = r * cos(th);
    yunit = r * sin(th);
    plot(xunit, yunit)
    plot(Cavity_points(:,1,i),Cavity_points(:,2,i))
    
    
    xlabel('X (m)')
    ylabel('Y (m)')
    legend('Initial Cavity','Deformed Cavity');
    
    figure(20)
   patch(fullcav(:,1),fullcav(:,2),'black')
   xlim([-1.4 1.4])
   ylim([-1.4 1.4])
   %axis off 

   F=getframe;
   %close all
   [X,~] = frame2im(F) ;
   XX=X(:,:,1)<1;
   stats{i} = regionprops(XX,'Eccentricity','MajorAxisLength','MinorAxisLength','Solidity');
   stats{i}.MinorAxisLength=stats{i}.MinorAxisLength*(1.4/(length(XX(:,2))));
   stats{i}.MajorAxisLength=stats{i}.MajorAxisLength*(1.4/(length(XX(:,1))));
    stats{i}.Area = polyarea(fullcav(:,1),fullcav(:,2));
    stats{i}.Alt_MajorAxis =C_Major_axis;
    stats{i}.Alt_MinorAxis =C_Minor_axis;
end

