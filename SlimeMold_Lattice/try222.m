%try normal power

m=10;
d=2;

k=1;

for i=1:100000
    z(i)=k*(norminv(rand(),m,d))^4;
end

histogram(z)