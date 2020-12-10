Xc = Grids{i,1}; Yc = Grids{i,2};


b = find(sign(u) ~= sign(Xc-50.1));

uc = u ; uc(b) = -1.*uc(b);
vc = v ; vc(b) = -1.*vc(b);

b = find(and(Yc>50.001,sign(vc)<0));

uc(b) = -1.*uc(b); vc(b) = -1.*vc(b);



ucc = inpaint_nans(uc);
vcc = inpaint_nans(vc);

figure;
subplot(1,2,1)
quiver(Xc,Yc,ucc,vcc)
axis image;
set(gca, 'YDir','reverse')

[sxv,syv] = meshgrid([-0.2 1]+50,1:2:70);
[sxh,syh] = meshgrid(0:2:100,[-0.1 1.1]+50);

subplot(1,2,2)

streamline(stream2(Xc,Yc,ucc,vcc,[sxv(:);sxh(:)],[syv(:);syh(:)]));
axis image;
set(gca,'YDir','reverse')




lx = linspace(4, 40, 40);
ly = linspace(60,50,40);

figure;
streamline(stream2(Xc,Yc,uc,vc,[lx;lx],[ly+2;ly-2]));
set(gca, 'YDir','reverse')

v0 = abs(v)<0.05;
u0 = abs(u)<0.05;

sx = [Xc(u0); Xc(v0)] ; sy = [Yc(u0); Yc(v0)];







            if iv==1
                b = find(and(Yc>50.001,sign(vc)<0));
                [sxv,syv] = meshgrid([-0.2 1]+50,1:2:70);
                [sxh,syh] = meshgrid(0:2:100,[-0.1 1.1]+50);
            else
                b = find(and(Yc>50.001,sign(vc)>0)); 
                [sxv,syv] = meshgrid([2 98],1:2:70);
                [sxh,syh] = meshgrid(0:2:100,[1 69]);
            end
            
            
            

