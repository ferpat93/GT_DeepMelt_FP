close all

Data =readmatrix('/Users/lfp3/Dropbox (GaTech)/GT/Spring-20/AutoAbaqus/Trial_Elastic_Disp_4.out','FileType','text');
% U_cav, X, Y, U1, U2, F1, F2
Data(:,8) = (Data(:,6).^2 + Data(:,7).^2).^(1/2);

Ps = unique(Data(:,1));
np = size(Data,1)/numel(Ps);
ang = linspace(0,180,np);

keys = 1:np:size(Data,1);

figure; hold on

for i=1:numel(keys)
    
    Data_i = sortrows(Data(keys(i):keys(i)+np-1,:),3);
    XY = Data_i(:,2:3);
    figure(1); hold on
    plot(XY(:,1),XY(:,2))
    
    figure(2); hold on
    plot(ang-90,Data_i(:,8))
    
end


[~,IND] = sortrows(Data(keys(i):keys(i)+np-1,:),2);
to_plot = linspace(1,np,5); IND = IND(to_plot);

for i=1:numel(IND)

    Data_i = Data(keys+IND(i)-1,[1 8]);

    figure(3); hold on
    plot(Data_i(:,1),Data_i(:,2))
    
end


figure;
P = abs(Data(:,1)); A = abs(Data(:,2)); S = abs(Data(:,3));
plot(A,P,'-o')
hold on

%plot([A(1) A(1,2)],[0 max(Data(:,1))])

d = 1800;
H = 30;
phi_d = 43.57; phi_r = deg2rad(phi_d);
ko = 1-sin(phi_r);
sv = 9.81*d*H/1000; sh = ko*sv;

plot([min(A) max(A)],[sv sv])
plot([min(A) max(A)],[sh sh])
xlabel('Cavity Area [m^2]')
ylabel('Pressure [kPa]')

figure; 
subplot(3,1,1)
plot(P,A); grid on
ylabel('Cavity Area [m^2]')
xlabel('Pressure [kPa]')

y=P; x=A;
yd = diff(y)./diff(x);
xd = (y(2:end)+y(1:(end-1)))/2;
subplot(3,1,2)
plot(xd,yd); grid on
%xlim([0 100])
xlabel('Pressure [kPa]')
ylabel('\delta P/ \delta A ')


subplot(3,1,3);
plot(P,S); grid on
xlabel('Pressure [kPa]')
ylabel('Aspect Ratio (S)')



