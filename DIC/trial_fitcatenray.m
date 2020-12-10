
%% fit catenary
close all

x = 0:0.25:10;

a = 7.5;

y = (a .* cosh((x-5)./a)) + (rand(size(x))-0.5)./200000000000000000;

scatter(x,y)


    Eqn = 'a*cosh((x-b)/a)+c';
    [f1,gn] = fit(x,y,Eqn);
    coeff = coeffvalues(f1);
    a = coeff(1); b = coeff(2); c = coeff(3);
    rsq = gn.rsquare;
