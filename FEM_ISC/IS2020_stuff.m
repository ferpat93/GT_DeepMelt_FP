%% Process Cavity Wall Data

% Load Struct containing Cavity points (16, one every 6 degrees

load('Cavity_Array.mat')
% Loads InitialCoord_Cav -> Initial coordinates of the cavity
% Loads Cav_coord -> 16X2X625 array with expanded cavity coordinates

load('CombinatorialCases.mat') 

load('EP_Array.mat') 

% E , v, phi, psi, Ko, G
CombinatorialCases = [CombinatorialCases 1-sin(deg2rad(CombinatorialCases(:,3))) CombinatorialCases(:,3)./(1+CombinatorialCases(:,2))];

% Order the points Counter-clock wise

InitialCoord_Cav = [InitialCoord_Cav; InitialCoord_Cav(2,:)];
InitialCoord_Cav(2,:) = [];

Cav_coord = [Cav_coord; Cav_coord(2,:,:)];
Cav_coord(2,:,:) = [];

figure()
plot(InitialCoord_Cav(:,2),InitialCoord_Cav(:,3))
hold on
for i=1:size(Cav_coord,3)
    plot(Cav_coord(:,1,i),Cav_coord(:,2,i))
end
axis equal


% Axis 
nT = size(Cav_coord,3);
Axis = zeros(nT,2); % Semi axis: Major/Minor
Areas = zeros(nT,2); % Poly, Ellipse
Perim = zeros(nT,2); % Poly, Ellipse

for i=1:nT
    Axis(i,:) = [Cav_coord(1,2,i) Cav_coord(end,1,i)];
    Areas(i,:) = [4*polyarea([Cav_coord(:,1,i);0],[Cav_coord(:,2,i);0]) pi*Axis(i,1)*Axis(i,2)];
    Perim(i,:) = Perimeters(Cav_coord(:,:,i),Axis(i,:));
end

AR = Axis(:,2)./Axis(:,1);

figure; scatter(Areas(:,1),Areas(:,2))
figure; scatter(Perim(:,1),Perim(:,2))

Area_Rsq = 1 - sum((Areas(:,1) - Areas(:,2)).^2)/sum((Areas(:,1) - mean(Areas(:,1))).^2);
Perim_Rsq = 1 - sum((Perim(:,1) - Perim(:,2)).^2)/sum((Perim(:,1) - mean(Perim(:,1))).^2);

% Indexes

% Extent
% 1. Area 2. A-Ao  3. MeanAxis 4. Eq.Diam 
IA = [Areas(:,1) Areas(:,1)-pi mean(Axis,2) (Areas(:,1)./pi).^(1/2)];

% Shape
% 1. AR 2. Circularity 3. Roundness 4.Axis diff 5. Perim
IS = [AR Areas(:,1)./(Perim(:,1).^2) Areas(:,1)./(Axis(:,1).^2) (Axis(:,1)-Axis(:,2)) Perim(:,1)];


K =zeros(nT,5);
for i=1:nT
    K(i,Y(i))=1;
end


X = [IA(:,1),IS(:,3)];
hist3(X,'CdataMode','auto')
xlabel('MPG')
ylabel('Weight')
colorbar
view(2)



%% Auxiliar Functions

function [P] = Perimeters(points,axis)

% Perimeter of polyshape
d = diff(points);
pp = 4*sum(sqrt(sum(d.*d,2)));

% Perimeter of Ellipse

pe = pi*((3*(sum(axis)))-((3*axis(1)+axis(2))*(axis(1)+3*axis(2)))^(1/2));

P = [pp pe];

end