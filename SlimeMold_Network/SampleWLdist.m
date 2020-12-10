

% Width
m = [-0.4 0.1 -0.2];
s = [0.61 0.61 0.55];
T = {'Control','Glucose','NaCl'};
w = linspace(0,2.5,50);

figure(1)
hold on
for t=1:3
    pd = makedist('Loglogistic','mu',m(t),'sigma',s(t));
    plot(w,pdf(pd,w));
end
legend(T)
xlabel('Width [mm]')
ylabel('pdf')



% Length
m = [0.37 0.6 0.22];
s = [0.49 0.49 0.55];
T = {'Control','Glucose','NaCl'};
w = linspace(0,6,50);

figure(2)
hold on
for t=1:3
    pd = makedist('Loglogistic','mu',m(t),'sigma',s(t));
    plot(w,pdf(pd,w));
end
legend(T)
xlabel('Length [mm]')
ylabel('pdf')