%% Cavity Analysis

% Cavity  
struct = load('EP_Array'); % Without Model 2! Erased bc didn't work
EP_Coord=struct.EP_Coord;
EP_Map=struct.EP_Map;

load('Cavity_Array')
load('Cavity_Stats')
a=cellfun(@(x) struct2cell(x),stats,'UniformOutput',false);
Cav_Indexes=cell2mat(cellfun(@(x) cell2mat(x'),a,'UniformOutput',false));
Cav_Indexes(2,:)=[];
% Cav_Indexes: +Axis -Axis Ecc Solidity Area Alt_+Axis Alt_-Axis


nModels=length(EP_Map(1,:));
EP_indexes=zeros(nModels,2);
%zeroT=1e-3;
r0=1;

ExtentMesh=40; % All EP Boundaries fall within the 40X40 neighborhood of the cavity

figure(1)
hold on
r=1;
th = 0:pi/100:pi/2;
xunit = r * cos(th);
yunit = r * sin(th);
plot(xunit, yunit)

xlabel('X (m)')
ylabel('Y (m)')
%legend('EP Boundary','Cavity Wall');
axis equal
xlim([0 2*ExtentMesh])
ylim([0 2*ExtentMesh])
    
[xMesh,yMesh] = meshgrid(0:.15:ExtentMesh, 0:.15:ExtentMesh);
statsEP=zeros(nModels,7); %[Eccentricity,Solidity,MinorAxis,MajorAxis,Area,MaxDist,AngletoMax]
contourEP=cell(nModels,1);

for i=1:15%nModels
    %close all
    [FullPastic, zeroT] = ThresHold(EP_Map(:,i));
    F_interp = scatteredInterpolant(EP_Coord(:,1),EP_Coord(:,2),EP_Map(:,i));
    F_interp.Method = 'natural';
    v = [zeroT,zeroT];
    
    figure(1)
    hold on
    EP_Contour=transpose(contour(xMesh,yMesh,F_interp(xMesh,yMesh),v));

    
% Plot to see 3d map of plastic strain 
%figure(2)
%surf(xMesh,yMesh,F_interp(xMesh,yMesh))
    
    maxX=max(EP_Contour(2:end,1));
    maxY=max(EP_Contour(2:end,2));
    contourEP{i}=[ EP_Contour(2:end,1) EP_Contour(2:end,2)];
    
%     fullcav=([[EP_Contour(2:end,1);-1.*flip(EP_Contour(2:end,1));-1.*EP_Contour(2:end,1);flip(EP_Contour(2:end,1))],...
%          [EP_Contour(2:end,2);flip(EP_Contour(2:end,2));-1.*EP_Contour(2:end,2);-1.*flip(EP_Contour(2:end,2))]
%                     ]) ;                      
%    figure(3)        
%    axis equal
%    patch(fullcav(:,1),fullcav(:,2),'black')
%    %xlim([-maxX maxX])
%    %ylim([-maxY maxY])
%    
%    %axis off 
% 
%    F=getframe;
%    clf
%    %close all
%    [X,~] = frame2im(F) ;
%    XX=X(:,:,1)<1;
%    a = regionprops(XX,'Eccentricity','Solidity');
%    a.MinorAxisLength=maxX;
%    a.MajorAxisLength=maxY;
%    %a.MinorAxisLength=statsEP{i}.MinorAxisLength*(maxX/(length(XX(:,2))));
%    %a.MajorAxisLength=statsEP{i}.MajorAxisLength*(maxY/(length(XX(:,1))));
%    a.Area = polyarea(fullcav(:,1),fullcav(:,2)); 
%    [a.MaxDist,I] = max((EP_Contour(2:end,1).^2+EP_Contour(2:end,2).^2).^0.5);
%    a.AngletoMax = atand(EP_Contour(I+1,2)/EP_Contour(I+1,1));
%    statsEP(i,:)=struct2array(a);
end


function [FullPlastic, T] = ThresHold(Map)
    %Function that computes the "zero" threshold
    minMap=min(Map);
    
    if minMap==0 %Not full Plasticity
        T=1e-2;
        FullPlastic = 0;
    else 
        T=prctile(Map,90);
        FullPlastic = 1;
    end

end

