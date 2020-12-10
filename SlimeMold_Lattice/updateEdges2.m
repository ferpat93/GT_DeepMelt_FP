function [newedges]=updateEdges2(L,edges,qEdges,Qin,me,de,nit)


%Column 1 of edges contains r
%Column 2 of edges contains alpha
%k=alpha*L/r^4 ----> Constant within the model

gamma=1.8; %Constant of tube growth
Qh=abs(Qin);
k=L*edges(1,2)/(edges(1,1)^4); %Sets constant value k=pi()/(8*n*L)
Dbar=k*(me^4+6*me^2*de^2+3*de^4);
%Io=median(abs(qEdges));
Io=Qh/nit;
a=(Io/Qh)^gamma

Dn=edges(:,2)./Dbar; %Auxiliar vector with values Dij
qm=(qEdges/Io).^gamma;
fq=((1+a).*qm)./(1+a.*qm);

ts=1/nit;
newedges(:,2)=Dbar.*max(Dn+ts*(fq-Dn),0); % Computes new Dij values
newedges(:,1)=((L/k).*newedges(:,2)).^(1/4); %Updates r from D

figure(1)
histogram(Dn)
figure(2)
histogram(qm)
figure(3)
histogram(fq)
figure(4)
histogram(fq-Dn)
figure(5)
histogram(Dn+ts*(fq-Dn))

newedges(:,2)/edges(:,2)