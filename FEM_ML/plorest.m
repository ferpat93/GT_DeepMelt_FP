close all

Data =readmatrix('/Users/lfp3/Dropbox (GaTech)/GT/Spring-20/AutoAbaqus/Pres_3/Pres_3.all','FileType','text');

Ps = unique(Data(:,1));
np = size(Data,1)/numel(Ps);

figure; hold on
%{
for i=1:numel(Ps)
    s = (i-1)*np+1;1)
    plot(,Data(s:s+np-1,3))
    
end
%}

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



