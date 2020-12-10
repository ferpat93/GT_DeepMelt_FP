%try normal power

m=0;
d=1;

k=1;

for i=1:1000
    z(i)=k*(norminv(rand(),m,d))^4;
end

histogram(z)