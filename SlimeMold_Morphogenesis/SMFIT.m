
function [SM, BLOB]=SMFIT(entities)

warning('off','all')
n=length(entities(:,1));
if n>500
    n=500;
end
x = linspace(1,n,n)';

Data=zeros(n,5,3);

for j=1:2:3

    y=entities(1:n,j);

    %Poly
    pdegree=0;
    pr2=0;
    while and(pr2<0.98,pdegree<9)
        pdegree=pdegree+1;
        fittype=strcat('poly',num2str(pdegree));
        [Pfitresult,gof] = fit(x,y,fittype); 
        pr2=gof.rsquare;
    end

    %Fourier
    fdegree=0;
    fr2=0;
    while and(fr2<0.95,fdegree<8)
        fdegree=fdegree+1;
        fittype=strcat('fourier',num2str(fdegree));
        [Ffitresult,gof] = fit(x,y,fittype); 
        fr2=gof.rsquare;
    end

    
    Py=feval(Pfitresult,x);
    Fy=feval(Ffitresult,x);
    Pdy=differentiate(Pfitresult, x);
    Fdy=differentiate(Ffitresult, x);

    Data(:,:,j)=[y Py Fy Pdy Fdy];
end

SM=Data(:,:,1);
BLOB=Data(:,:,3);

