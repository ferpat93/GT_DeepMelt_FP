%% Cavity Analysis

% Cavity  
EP_Folder='C:\Users\lfp3\Documents\IS_SensitivityAnalysis\Cases\EP\';
struct = load(fullfile(EP_Folder,'EP_Array'));
EP_Coord=struct.EP_Coord;
EP_Map=struct.EP_Map;

nModels=length(EP_Map(1,:));
EP_indexes=zeros(nModels,2);
zeroT=1e-2;
r0=1;

ExtentMesh=40;
[xMesh,yMesh] = meshgrid(0:.15:ExtentMesh, 0:.15:ExtentMesh);

for i=200:200%nModels           
    F_interp = scatteredInterpolant(EP_Coord(:,1),EP_Coord(:,2),EP_Map(:,i));
    F_interp.Method = 'natural';
    v = [zeroT,zeroT];
    EP_Contour=contour(xMesh,yMesh,F_interp(xMesh,yMesh),v);

    % Fit Cavity Ellipse and get axis
    ellipse_struct=fit_ellipse(EP_Contour(1,:),EP_Contour(2,:)); % Struct containing ellipse parameters
    EP_Major_axis=ellipse_struct.long_radius;
    EP_Minor_axis=ellipse_struct.short_radius;

    %Find Indexes
    EP_indexes(i,1)=1*(((EP_Minor_axis*EP_Major_axis)-(EP_Major_axis*EP_Minor_axis))/r0^2);
    EP_indexes(i,2)=EP_Minor_axis/EP_Major_axis; 
end
