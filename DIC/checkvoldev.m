
clc

x = 50 ; y = 15 ; z = 7 ;

u = reshape(permute(S(x,y,z,3:11),[4 3 2 1]),[3 3]);
up = u + eye(3);

vol = det(up)-1
volr = S(x,y,z,2)

Udev = up*(det(up)^(-1/3));

dev = norm(Udev-eye(3),2)
devr = S(x,y,z,1)

dev/devr