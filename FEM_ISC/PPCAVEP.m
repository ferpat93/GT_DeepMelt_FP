


% Fit Cavity Ellipse and get axis
ellipse_struct=fit_ellipse(Cavity_points(:,1),Cavity_points(:,2)); % Struct containing ellipse parameters
assignin('base','CavityPoints',Cavity_points)
assignin('base','Cavity_struct',ellipse_struct)

C_Major_axis=ellipse_struct.long_radius
C_Minor_axis=ellipse_struct.short_radius

%Fit EP Boundary Ellipse and get axis
ellipse_struct=fit_ellipse(Ellipse_points(:,1),Ellipse_points(:,2)); % Struct containing ellipse parameters

assignin('base','Ellipse_points',Ellipse_points)
assignin('base','EP_struct',ellipse_struct)

EP_Major_axis=ellipse_struct.long_radius
EP_Minor_axis=ellipse_struct.short_radius

%Find Indexes

%Cavity
Index(1,1)=100*(((C_Major_axis*C_Minor_axis)/(r0^2))-1);
Index(1,2)=C_Minor_axis/C_Major_axis;

%EP Boundary
Index(1,3)=1*(((EP_Minor_axis*EP_Major_axis)-(C_Major_axis*C_Minor_axis))/r0^2);
Index(1,4)=EP_Minor_axis/EP_Major_axis;