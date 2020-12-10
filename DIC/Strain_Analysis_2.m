
close all

%% Lines to use tight subplots
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.065 0.065], [0.2 0.05], [0.075 0.04]); % subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
if ~make_it_tight,  clear subplot;  end

load('Data1.mat');

TestVariables = [ [2 4 6 2 4 6 2 4 2 6 6 4 6]' [1 1 1 1 1 1 0 0 0 0 0 0 0]' [0 0 0 1 1 1 0 0 1 0 1 1 2]']; % Length (2,4,6), density (0,1), load(0,1,2)
Order_table = [7 8 10 9 12 11 13 ; 1 2 3 4 5 6 nan];

%% Principal strain orientation

% Get angles
getAngles = false;

iv = 1; % Major principal strain 
sy = 1;
components = [ 1 2 ; 3 2 ; 3 1];
sizes = [100 70 100];
    
if getAngles 
    angles = cell([size(Strain) 2]);
    for test = 1:size(angles,1) % Loop tests
        for step=1:size(angles,3) % Loop steps
            
            A = Strain{test,step}(:,:,:,12:end);
            ix = 100/size(A,1); xc = ix/2:ix:100-ix/2;
            yc = xc; iy = 70/40;%size(A,3);
            zc = 1.*(iy/2:iy:70-ix/2);

            [Grids{1,1}, Grids{1,2}] = meshgrid(xc,zc); % Front
            [Grids{2,1}, Grids{2,2}] = meshgrid(yc,zc); % Side
            [Grids{3,1}, Grids{3,2}] = meshgrid(xc,yc); % Top

            P = cat(4,A(:,:,:,iv*3:iv*3+2),abs(A(:,:,:,12)-A(:,:,:,14))); % Principal Strain Array      

            Slices{1} = flipud(permute(P(:,round(size(P,2)/2),1:40,:),[3 1 4 2])); % Front
            Slices{2} = flipud(permute(P(round(size(P,2)/2),:,1:40,:),[3 2 4 1])); % Side  
            Slices{3} = P(:,:,round(size(P,3)/2),:); % Top

            for view=1:3 % Front, Side and Top
                u = Slices{view}(:,:,components(view,1)); v = Slices{view}(:,:,components(view,2));
                shear = imresize(Slices{view}(:,:,4),[sizes(components(view,2)) sizes(components(view,1))]);

                Xc = Grids{view,1}; Yc = Grids{view,2};
                ys = sign(Yc-50.1); xs = sign(Xc-50.1); 

                b = (sign(u)==xs)+2.*(sy.*sign(v)==ys);
                uc = u ; vc = v ;
                uc(b==0) = -1.*uc(b==0); vc(b==0) = -1.*vc(b==0); % Change both (opposite)
                uc(or(b==1,b==2)) = nan; vc(or(b==1,b==2)) = nan; % Change trouble-some (nan)
                ucc = inpaint_nans(uc); vcc = inpaint_nans(vc); % Interpolate vector components

                angles{test,view,step} = rad2deg(abs(atan(vcc./ucc))); %rad2deg(atan2(vcc,ucc));
                
%                 figure; hold on
%                 A = rad2deg(atan2(vcc,ucc));
%                 rad2deg(atan(vcc./ucc))
%                 imagesc(A)
%                 quiver(Xc,Yc,ucc,vcc)
%                 gs=1;
%                 
%                 
                
            end
        
        end   
    end
    save('angles.mat','angles');    
else
    load('angles.mat')
end

Tests = [1 13]; step = 2;
%PlotPrincStrainOrientation(Strain,Tests,step)

%% Shear Bands
getBands = false;

if getBands
    
    %Bands = cell(size(Strain));
    lim = [0.05 0.1]; 
        
    for t=1:size(Strain,1)    
        for step=2:2
            %I = Porosity{t,3}-Porosity{t,1}; % Porosity Change
            I = Strain{t,step}(:,:,:,end); % Max Shear
            [~,Front,~] = get_slices(I);  
            Front(isnan(Front))=max(Front(:));
            Front = imresize(Front,round(1.75*size(Front)),'bilinear');
            %Front = imresize(Front,[84 100]); 
            
            %figure;
            figure('units','normalized','outerposition',[0 0 1 1]);
                imagesc(Front,lim);
                colormap(hot)
                colorbar
                
                %imcontour(Front_D,[0.02 0.1])
                %[C,h]=imcontour(Front_D,[0.02 0.1]);
                %clabel(C,h)
                %text(20,0,num2str(t))
                %F = getframe(gcf);
                %[X, Map] = frame2im(F);
                %imshow(X)
                axis image;
                h = drawfreehand;
                j = drawfreehand;
                Bands{t,step} = cat(3,createMask(h),createMask(j));
            close;

            save('bands.mat','Bands');
        end
    end 
    
else
    load('bands.mat');
end


%% Catenary
getCatenary = false;

if getCatenary
    
    %Catenary = nan(size(PS,1),4);
    lim = [0.025 0.1]; 
        
    for t=1:size(Strain,1)    
        for step=2:2

            [~,Front_D,~] = get_slices(Strain{t,step}(:,:,:,end));  % Shear
            Front_D(isnan(Front_D))=1;
            Front_D=Front_D(2:end,1:end-1);
            if size(Front_D)>40; rows = 84; else; rows = 70; end
            Front_D = imresize(Front_D,[rows 100]);
            
            %figure;
            figure('units','normalized','outerposition',[0 0 1 1]);
                imagesc(Front_D,lim);
                colormap(hot)
                axis image;
                contour = drawpolyline;
                x=contour.Position(:,1); y=90-contour.Position(:,2);
                [a,b,c,rsq] = fit_catenary(x,y);
                Catenary(t,:) = [ a b c rsq];
            close;
                
        end
    end 
    
    save('catenary.mat','Catenary');
else
    load('catenary.mat');
end

%% Analyze Catenary

catenary_table = nan([size(Order_table) 4]);
for t=1:size(Bands,1)     
    [r,c]=find(Order_table==t);
    catenary_table(r,c,:)=Catenary(t,:);
end

%% Analyze bands 
step = 2;
d50 = 0.35; %mm

wa_bands = nan([size(Order_table) 5]);

for t=1:size(Bands,1)  
    
    [r,c]=find(Order_table==t);
    band = Bands{t,step};

    %figure;
    for bs = 1:size(band,3)
        bi = band(:,:,bs);
        skel = bwskel(bi);
        
         %subplot(1,2,bs)
         %imagesc(bi+skel)
        
        if sum(skel(:))>10
            dist = bwdist(~bi);
            veins = skel.*dist;
            ws = veins(:); ws(ws==0)=[];
            wa_bands(r,c,1) = nanmean([2*median(ws) wa_bands(r,c,1)]);
            wa_bands(r,c,2) = wa_bands(r,c,1)/d50;
            
            [y,x] = find(skel); % Or skel
            X = [ones(length(x),1) x];
            b = X\y;
            
            [~,~,R,~] = circfit(x,y);
            wa_bands(r,c,5)=1/R;
            
            wa_bands(r,c,3) = nanmean([abs(atand(b(2))) wa_bands(r,c,3)]);
            
            xy = sortrows([x y],2,'descend');
            xy = xy(1:round(2*size(xy,1)/3),:);
            
            X = [ones(size(xy,1),1) xy(:,1)];
            b = X\xy(:,2);

            wa_bands(r,c,4) = nanmean([abs(atand(b(2))) wa_bands(r,c,4)]);
        
            fgddg=1;
        else
            Bands{t,step}(:,:,bs) = zeros(size(skel));
        end
        
    end
    


end

%% Get change inside and outside shear bands -> Porosity , deviatoric and Volumetric Strain

load('NewPorosity.mat');

porosity_volumetric_table = nan([size(Order_table) 6]); % Deviatoric (O,I) Volumetric (O,I) Porosity(O,I)

for t=1:13
    
        [r,c]=find(Order_table==t); 
        
        % Porosity change: Initial vs. 2nd inlfation
        
        [~,Front_P,~] = get_slices(Porosity{t,3}-Porosity{t,1});  % Get Front and Side View
        Front_P = imresize(Front_P,[84 100]); 
        
        [~,Front_D,~] = get_slices(Strain{t,2}(:,:,:,1));  % Deviatoric
        Front_D(isinf(Front_D))=nan; Front_D = inpaint_nans(Front_D);
        Front_D = imresize(Front_D,round(1.75*size(Front_D)));
        
        [~,Front_V,~] = get_slices(Strain{t,2}(:,:,:,2));  % Volumetric   
        Front_V(isinf(Front_V))=nan; Front_V = inpaint_nans(Front_V);        
        Front_V = imresize(Front_V,round(1.75*size(Front_V)));
        % Store porosity inside and outside
        band = sum(Bands{t,2},3); 
        
        % Deviatoric
        porosity_volumetric_table(r,c,1)=100*nanmean(Front_D(~logical(band))); % Outside  
        porosity_volumetric_table(r,c,2)=100*nanmean(Front_D(logical(band))); % Band  
      
        % Volumetric
        porosity_volumetric_table(r,c,3)=100*nanmean(Front_V(~logical(band))); % Outside  
        porosity_volumetric_table(r,c,4)=100*nanmean(Front_V(logical(band))); % Band  
      
        % Porosity      
        band(:,101:end)=[];
        if size(band,1)<size(Front_P,1); band = [band ; zeros(14,100)]; end
        porosity_volumetric_table(r,c,5)=nanmean(Front_P(~logical(band))); % Outside  
        porosity_volumetric_table(r,c,6)=nanmean(Front_P(logical(band))); % Band  

        sgss=1;
        
end

figure; hold on
scatter(porosity_volumetric_table(1,:,4),porosity_volumetric_table(1,:,6),'r')
scatter(porosity_volumetric_table(2,:,4),porosity_volumetric_table(2,:,6),'b')

%% Terminal Void ratio

VoidRatios = zeros(13,3);
Terminal_e_table = zeros(size(Order_table));

for t=1:13
    band = sum(Bands{t,2},3); 
    band(:,101:end)=[];
        
    for step=1:3
        
        [~,Front,~] = get_slices(Porosity{t,step});  % Get Front View  
        Front = imresize(Front,[84 100]); 
        e = Front ./ (100-Front);
        
        if size(band,1)<size(Front,1); band = [band ; zeros(14,100)]; end
        
        if step==1
            VoidRatios(t,step)=nanmean(e(:)); % Band   
        else
            VoidRatios(t,step)=nanmean(e(logical(band))); % Band      
        end
    end
    
    [r,c]=find(Order_table==t); 
    Terminal_e_table(r,c)=VoidRatios(t,step);
end

VoidRatios(13,3)=VoidRatios(13,2)-0.0542;

for t=1:13
    [r,c]=find(Order_table==t); 
    Terminal_e_table(r,c)=VoidRatios(t,step);
end

colorlines = {'b','k','r'};
linetyp = {'-','-.',':'};

figure;
% Loose
subplot(1,2,1); hold on
ts = [8 10 12 11 13];
for ti=1:numel(ts)
    t = ts(ti);
    plot(1:3,VoidRatios(t,1:3),[colorlines{TestVariables(t,1)/2} linetyp{TestVariables(t,3)+1}])           
end

legend({'L = 4 - \sigma1','L = 6 - \sigma1','L = 4 - \sigma2','L = 6 - \sigma2','L = 6 - \sigma3'},'Location','southoutside','NumColumns',3);
xlim([0.8 3.2])
ylabel('Void ratio (e)')
title('Loose Specimens')

subplot(1,2,2); hold on
ts = [1 2 3 4 6];
for ti=1:numel(ts)
    t = ts(ti);
    plot(1:3,VoidRatios(t,1:3),[colorlines{TestVariables(t,1)/2} linetyp{TestVariables(t,3)+1}])   
end

legend({'L = 2 - \sigma1','L = 4 - \sigma1','L = 6 - \sigma1','L = 2 - \sigma2','L = 6 - \sigma2'},'Location','southoutside','NumColumns',2);
xlim([0.8 3.2])
ylabel('Void ratio (e)')
title('Dense Specimens')

%% PLOT - Porosity

divMap = [brewermap(101,'RdBu')]; %Divergent Colormap (- to +)
seqMap = [brewermap(101,'YlGn')]; %Sequential Colormap (0 to inf)

testsToPlot = [9 12 11; 4 5 6]; % Tests to plot rows: loose/dense

maplim = [15 25];
titles = {'L = 2D' , 'L = 4D' , 'L = 6D' };
for di=1:size(testsToPlot,1) % Loop over density
    
    figure;
    
    for li = 1:size(testsToPlot,2) % Loop over lengths
        P = cell(3,1);
        t = testsToPlot(di,li);

        P{1} = Porosity{t,1}; % Initial Porosity
        P1 = Porosity{t,2}; P2 = Porosity{t,3};
        P{2} = P1-P{1}; % Porosity change: Initial vs. 1st inflation
        P{3} = P2-P{1}; % Porosity change: Initial vs. 2nd inlfation

        [~,Front,Side] = get_slices(P{3});  % Get Front and Side View
        Front = imresize(Front,[84 100]); Side = imresize(Side,[84 100]);
        
        %maglim = max(abs(prctile([Front(:);Side(:)],5)),prctile([Front(:);Side(:)],95));
        maglim = maplim(di);
        lim = [-maglim maglim];
        ax(1)=subplot(2,3,li); imagesc(fliplr(Front),lim); colormap(ax(1),divMap);axis image ; axis off;
        title(titles{li});
        ax(1+3)=subplot(2,3,li+3); imagesc(fliplr(Side),lim); colormap(ax(1+3),divMap);axis image ; axis off;

        % Old Code to show initial porosity and 2 steps for each test
%         [~,Front,Side] = get_slices(P{i});  % Get Front and Side View
%         %maglim = max(abs(extremes(t,2,step)),extremes(t,3,step)); lim = [-maglim maglim];
%         maglim = max(abs(prctile([Front(:);Side(:)],5)),prctile([Front(:);Side(:)],95))
%         if i==1; lim = [0 100]; CM =  seqMap; else; lim = [-maglim maglim]; CM =  divMap;  end
%         ax(i)=subplot(2,3,i); imagesc(Front,lim); colormap(ax(i),CM);axis image ; axis off;
%         ax(i+3)=subplot(2,3,i+3); imagesc(Side,lim); colormap(ax(i+3),CM);axis image ; axis off;
%         colorbar('southoutside')
    
    end
    
    c=colorbar('east'); c.Ruler.TickLabelFormat='%g%%';
    h = text(0.1,0.5,'Front'); set(h,'Rotation',90);
    h = text(0.1,0.5,'Side'); set(h,'Rotation',90);

end 

%% PLOT - Volumetric vs Deviatoric Strain
%%% -> Compressive volumetric strains are negative
% NOT USED FROM NOW BC EXTREMES ARE NOT IN SECTIONS
%A = cellfun(@(x) getmaxmin(x),Strain(:,1:2),'UniformOutput',false);
%extremes(:,:,1)=cell2mat(A(:,1));
%extremes(:,:,2)=cell2mat(A(:,2)); % Max Deviatoric - Min Volumetric - Max Volumetric

% Plot to check strain maps

step=2;
divMap=[[0 0 0]; brewermap(101,'RdBu')];
seqMap=[brewermap(101,'YlGn');[0 0 0]];

testsToPlot = [9 12 11; 4 5 6]; % Tests to plot rows: loose/dense

devlims = [10 10]; vollims = [2.5 5];

for di=1:size(testsToPlot,1) % Loop over density
    
    Dev = figure;
    Vol = figure;
    
    for li = 1:size(testsToPlot,2) % Loop over lengths
    
    t = testsToPlot(di,li);

    %Deviatoric
    figure(Dev)
    [~,Front_DD,Side_DD] = get_slices(Strain{t,step}(:,:,:,1));  % Deviatoric
    Front_D = 100.*Front_DD(1:40,:) ; Side_D = 100.*Side_DD(1:40,:);
    
    %lim = [0 extremes(t,1,step)*1.01];
    lim = [ 0 1.01*prctile([Front_D(isfinite(Front_D));Side_D(isfinite(Side_D))],97.5)];
    lim = [0 devlims(di)];
    ax(1)=subplot(2,3,li); imagesc(Front_D,lim); colormap(ax(1),seqMap);axis image ; axis off;
        ax(4)=subplot(2,3,li+3); imagesc(Side_D,lim); colormap(ax(4),seqMap);axis image ; axis off;
    colorbar('southoutside')
    
    % Volumetric
    figure(Vol)
    [~,Front_VV,Side_VV] = get_slices(Strain{t,step}(:,:,:,2));  % Volumetric
    Front_V = 100.*Front_VV(1:40,:) ; Side_V = 100.*Side_VV(1:40,:);
    %maglim = max(abs(extremes(t,2,step)),extremes(t,3,step)); lim = [-maglim maglim];
    maglim = max(abs(prctile([Front_V(isfinite(Front_V));Side_V(isfinite(Side_V))],5)),prctile([Front_V(isfinite(Front_V));Side_V(isfinite(Side_V))],95));
    lim = [-maglim maglim];
    lim = [-vollims(di) vollims(di)];
    ax(2)=subplot(2,3,li); imagesc(Front_V,lim); colormap(ax(2),divMap);axis image ; axis off;
    ax(5)=subplot(2,3,li+3); imagesc(Side_V,lim); colormap(ax(5),divMap);axis image ; axis off;
    colorbar('southoutside')
    
    % Ratio
%     FR = -1*Front_D./Front_V;
%     SR = -1*Side_D./Side_V;
%     
%     maglim = max(abs(prctile(FR(:),10)), abs(prctile(FR(:),90))); lim = [-maglim maglim];
%     ax(3)=subplot(2,3,3); imagesc(FR,lim); colormap(ax(3),divMap);axis image ; axis off;
%     maglim = max(abs(prctile(SR(:),10)), abs(prctile(SR(:),90))); lim = [-maglim maglim];
%     ax(6)=subplot(2,3,6); imagesc(SR,lim); colormap(ax(6),divMap);axis image ; axis off;
%     colorbar('southoutside')
        
    %[~,Front_M,Side_M] = get_slices(MP{t,step}(:,:,:,1));  % Magnitude
    
    end
end 

%% Plain Strain Soil
% (Probe axis is along Y - 2nd dimension)
% Norm difference (2D vs. 3D sensor), figures alongside percentages of
% variation.

v=0:(100/58):100;
pPS = zeros(13,2);
testsToPlot = [9 12 11];
ImagesToPlot = cell(2,numel(testsToPlot),2);
maxStrain = [0 0]; steptoplot = 1;

colorlines = {'r','k','b'};
linetyp = {'-.','-'};

for l=2:2:6
    
    idx = find(TestVariables(:,1)==l);
    intv = (-l*5)-5:1.25:(l*5)+5;
    
    for i=1:numel(idx)

        D = Device{idx(i),1}; center = 0.07*(D(1,1)+D(end,1))/2;
        if abs(center-50)>5; center =50; end
        vc = v - center;
            
        for s=1:2
            
            Loss = PS{idx(i),s}(:,:,:,2); % Strain Loss 
            Mag = PS{idx(i),s}(:,:,:,1); % Strain Magnitude
            %meanLoss = nansum(Loss,[1,3])./sum(~isnan(Loss),[1,3]);
            meanLoss = nansum(Loss.*Mag,[1,3])./nansum(Mag,[1,3]);
            intLoss = interp1(vc,meanLoss,intv); 
            
            figure(s+35)
            subplot(1,3,l/2)
            hold on
            plot(intv,intLoss,[colorlines{TestVariables(idx(i),3)+1} linetyp{TestVariables(idx(i),2)+1}])
            
            pPS(idx(i),s) = 100*(sum(abs(intLoss)<=5)/sum(abs(intv)<=l*5));
            
        end
        
        ip = find(testsToPlot==idx(i));
        if ~isempty(ip)   
            psi = Loss; 
            maxStrain(1) = max(maxStrain(1),max(Mag(~isinf(Mag)))); maxStrain(2) = max(maxStrain(2),max(psi(~isinf(psi))));
            Mag(isnan(Mag)) = inf;   psi(isnan(psi)) = inf;       
            [~,ImagesToPlot{1,ip,1},ImagesToPlot{2,ip,1}] = get_slices(Mag);
            [ImagesToPlot{1,ip,2},~,ImagesToPlot{2,ip,2}] = get_slices(psi); 
            
        end

    end
    
    figure(35+1); subplot(1,3,l/2); xlim([min(intv) max(intv)]); ylim([0 60]); xlabel('Distance from axis center [mm]'); ylabel('Average Plane Strain Deviation [%]'); title(strcat('L =',{' '},num2str(l),'D')); grid on
    figure(35+2); subplot(1,3,l/2); xlim([min(intv) max(intv)]); ylim([0 45]); xlabel('Distance from axis center [mm]'); ylabel('Average Plane Strain Deviation [%]'); title(strcat('L =',{' '},num2str(l),'D')); grid on
end

for s=1:2; figure(s+35); subplot(1,3,3); legend({'Dense - s1','Dense - s2','Loose - s1','Loose - s2','Loose - s3'},'Location','southoutside'); end

% Plot Strain Magnitude and Loss

divMap=[brewermap(100,'Oranges');[0 0 0]];
maxStrain = [0.1 100];

for li=1:3
    figure(20+1)
    subplot(2,3,li); imagesc(ImagesToPlot{1,li,1}(1:40,:),[0 maxStrain(1)]);  axis image ; axis off;
    title(strcat('L = ',num2str(2*li),'D'))
    subplot(2,3,li+3); imagesc(ImagesToPlot{2,li,1}(1:40,:),[0 maxStrain(1)]); axis image ; axis off;
    colormap(divMap)  
    
    figure(20+2)
    subplot(2,3,li); imagesc(ImagesToPlot{1,li,2},[0 maxStrain(2)]);  axis image ; axis off;
    title(strcat('L = ',num2str(2*li),'D'))
    subplot(2,3,li+3); imagesc(ImagesToPlot{2,li,2},[0 maxStrain(2)]); axis image ; axis off;
    colormap(divMap)  
    
end 

figure(20+1)
    colorbar('eastoutside')
    h(1)=text(0.05,0.7,'Transversal Section'); h(2)=text(-2.65,0.2,'Longitudinal Section');
    set(h(1),'Rotation',90); set(h(2),'Rotation',90);
figure(20+2)
    colorbar('eastoutside')
    h(1)=text(0.05,0.7,'Horizontal Section'); h(2)=text(-2.65,0.2,'Longitudinal Section');
    set(h(1),'Rotation',90); set(h(2),'Rotation',90);

% Difference between inflation steps

dif = pPS(:,2)-pPS(:,1);
xv = (TestVariables(:,1))/2 + (TestVariables(:,2)/2) ;

figure; h=gscatter(xv,dif,TestVariables(:,3),'rkb','od*');
ax = gca; ax.XTick = unique(sort(xv));
DensityLabels = {'Loose','Dense','Loose','Dense','Loose','Dense'}';
ax.XTickLabels = DensityLabels;
grid on; ylabel({'Change of length in'; 'plane Strain [%]'})
text(1.17,36,'L = 2D'); text(2.17,36,'L = 4D'); text(3.17,36,'L = 6D')
hold on; k = plot([1.75 1.75 2.75 2.75],[-20 100 100 -20]);
ylim([-10 40]); xlim([0.75 3.75]); legend(h,'s1','s2','s3')


% Plot the percentages by variable

xv = (TestVariables(:,1))/2 + (TestVariables(:,2)/2) ;

figure; h=gscatter(xv,pPS(:,2),TestVariables(:,3),'rkb','od*');
ax = gca; ax.XTick = unique(sort(xv));
DensityLabels = {'Loose','Dense','Loose','Dense','Loose','Dense'}';
ax.XTickLabels = DensityLabels;
grid on; ylabel({'Length of Cavity in';' Plane Strain [%]'})
text(1.17,82,'L = 2D'); text(2.17,82,'L = 4D'); text(3.17,82,'L = 6D')
hold on; k = plot([1.75 1.75 2.75 2.75],[20 100 100 20]);
ylim([40 85]); xlim([0.75 3.75]); legend(h,'s1','s2','s3')

    
%% OLD STUFF %%

for t=1:numel(tests)
    test = tests(t);
    
    S = Results{test};
    
    step = 2;
    Disp_step2=S.Disp{step}{4};
    Disp_step2(Disp_step2==Inf)=-Inf;
    max_mag = 0.07*max(Disp_step2(:));
    
    [~,Front,~] = get_slices(S.Disp{step}{3});
    abs_front = 0.07*abs(Front);
    abs_front(Front==-Inf)=nan;
    
    figure(40+t)
    subplot(1,2,1)
    imagesc(abs_front)
    axis image
    
    caxis([0 max_mag])
    Cmap= [ flipud(hot(50)) ; 0 0 1];
    title('Y Displacement [mm]')
    colormap(Cmap)
    %colorbar
    
    [~,Front,~] = get_slices(S.Disp{step}{4});
    
    figure(40+t)
    subplot(1,2,2)   
    imagesc(0.07.*Front)
    title('Displacement (magnitude) [mm]')
    axis image
    colorbar
    
end

%% Effect of variables

% Length

tests = [];
%tests = [1 2 3];

length_d = {'L = 2 cm','L = 4 cm','L = 6 cm'};

for t=1:numel(tests)
    test = tests(t);    
    S = Results{test};
    Device = S.Device;
    
    length_array = size(Device{1},2);
    crop_length = 2.5; % On both sides
    vp = [round(length_array*(crop_length/100)) ; length_array - round(length_array*(crop_length/100))];
    xv = linspace(crop_length,100-crop_length,vp(2)-vp(1)+1);       
    
    for s=1:3%numel(Device)
        Indexes = get_probe_geometry(Device{s});
        
        figure(10+t)
        subplot(1,3,1)
        hold on
        plot(xv,Indexes(1,vp(1):vp(2)))
        title('Cross section area')
        xlim([0 100])
        ylim([80 180])
        xlabel('Percentage of length') 
        ylabel('Area [mm2]') 
        
        figure(10+t)
        subplot(1,3,2)
        hold on
        plot(xv,Indexes(2,vp(1):vp(2)))
        title('Major Axis length')
        xlim([0 100])
        ylim([10 16])
        xlabel('Percentage of length') 
        ylabel('Length [mm]') 
        
        figure(10+t)
        subplot(1,3,3)
        hold on
        plot(xv,Indexes(2,vp(1):vp(2))./Indexes(3,vp(1):vp(2)));
        title('Eccentricity')
        xlim([0 100])
        ylim([1 1.15])
        xlabel('Percentage of length') 
        ylabel('Ratio of major/minor axis') 
        
        if s==3
            figure(100)
            subplot(1,3,1)
            hold on
            plot(xv,Indexes(1,vp(1):vp(2)))
            title('Cross section area')
            xlim([0 100])
            ylim([80 180])
            xlabel('Percentage of length') 
            ylabel('Area [mm2]') 

            subplot(1,3,2)
            hold on
            plot(xv,Indexes(2,vp(1):vp(2)))
            title('Major Axis length')
            xlim([0 100])
            ylim([10 16])
            xlabel('Percentage of length') 
            ylabel('Length [mm]') 

            subplot(1,3,3)
            hold on
            plot(xv,Indexes(2,vp(1):vp(2))./Indexes(3,vp(1):vp(2)));
            title('Eccentricity')
            xlim([0 100])
            ylim([1 1.15])
            xlabel('Percentage of length') 
            ylabel('Ratio of major/minor axis') 
        end 
    end
    
    figure(100)
    legend(length_d)
    
    figure(10+t)
    legend('P = 0 kPa','P = 100 kPa','P = 200 kPa')
    sgtitle(['Effect of length - ' length_d{t}])
end


% Density

tests = [];
%tests = [2 8 3 10];

% Per test
for t=1:numel(tests)
    test = tests(t);    
    S = Results{test};
    Device = S.Device;
    
    length_array = size(Device{1},2);
    crop_length = 2.5; % On both sides
    vp = [round(length_array*(crop_length/100)) ; length_array - round(length_array*(crop_length/100))];
    xv = linspace(crop_length,100-crop_length,vp(2)-vp(1)+1);       
    
    for s=1:numel(Device)
        Indexes = get_probe_geometry(Device{s});
        
        figure(20+t)
        subplot(1,3,1)
        hold on
        plot(xv,Indexes(1,vp(1):vp(2)))
        title('Cross section area')
        xlabel('Percentage of length') 
        ylabel('Area [mm2]') 
        
        figure(20+t)
        subplot(1,3,2)
        hold on
        plot(xv,Indexes(2,vp(1):vp(2)))
        title('Major Axis')
        xlabel('Percentage of length') 
        ylabel('Axis length [mm]') 
        
        figure(20+t)
        subplot(1,3,3)
        hold on
        plot(xv,Indexes(2,vp(1):vp(2))./Indexes(3,vp(1):vp(2)));
        title('Ratio between Axis')
        xlabel('Percentage of length') 
        
    end
    
    legend('Step 1','Step 2','Step 3')
    sgtitle(['Effect of Density - ' tests_list{test}])
    
end

tests = [ 7 1; 8 2; 10 3];
tests = [];

% Summary of tests
for t=1:size(tests,1)
    
    for c=1:2
    test = tests(t,c);    
    S = Results{test};
    Device = S.Device;
    
    length_array = size(Device{1},2);
    crop_length = 2.5; % On both sides
    vp = [round(length_array*(crop_length/100)) ; length_array - round(length_array*(crop_length/100))];
    xv = linspace(crop_length,100-crop_length,vp(2)-vp(1)+1);       
    
    s = 3;
      
        
    Indexes = get_probe_geometry(Device{s});
    
        figure(101)
        subplot(1,3,t)
        hold on
        plot(xv,Indexes(1,vp(1):vp(2)))
        title(length_d{t})
        xlim([0 100])
        %ylim([80 200])
        xlabel('Percentage of length') 
        ylabel('Area [mm2]') 

        figure(102)
        subplot(1,3,t)
        hold on
        plot(xv,Indexes(2,vp(1):vp(2)))
        title(length_d{t})
        xlim([0 100])
        %ylim([10 20])
        xlabel('Percentage of length') 
        ylabel('Length [mm]') 

        figure(103)
        subplot(1,3,t)
        hold on
        plot(xv,Indexes(2,vp(1):vp(2))./Indexes(3,vp(1):vp(2)));
        title(length_d{t})
        xlim([0 100])
        %ylim([1 1.2])
        xlabel('Percentage of length') 
        ylabel('Ratio of major/minor axis') 
        
    end
end   

sig={'Dr = 40%','Dr = 75%'};
figure(101)
legend(sig)
sgtitle('Cross section area [mm2]')
figure(102)
legend(sig)
sgtitle('Major Axis length [mm]')
figure(103)
legend(sig)
sgtitle('Eccentricity')

% Single test L6
tests = [10 3];
tests = [];

for c=1:numel(tests)
    test = tests(c);    
    S = Results{test};
    Device = S.Device;
    
    length_array = size(Device{1},2);
    crop_length = 2.5; % On both sides
    vp = [round(length_array*(crop_length/100)) ; length_array - round(length_array*(crop_length/100))];
    xv = linspace(crop_length,100-crop_length,vp(2)-vp(1)+1);       
    
    s = 3;
         
    Indexes = get_probe_geometry(Device{s});
    
    figure(111)
    subplot(1,3,1)
    hold on
    plot(xv,Indexes(1,vp(1):vp(2)))
    xlim([0 100])
    %ylim([80 200])
    xlabel('Percentage of length') 
    ylabel('Area [mm2]') 
    title('Cross section area [mm2]')

    subplot(1,3,2)
    hold on
    plot(xv,Indexes(2,vp(1):vp(2)))
    xlim([0 100])
    %ylim([10 20])
    xlabel('Percentage of length') 
    ylabel('Length [mm]') 
    title('Major Axis length [mm]')
    
    subplot(1,3,3)
    hold on
    plot(xv,Indexes(2,vp(1):vp(2))./Indexes(3,vp(1):vp(2)));
    xlim([0 100])
    %ylim([1 1.2])
    xlabel('Percentage of length') 
    ylabel('Ratio of major/minor axis') 
    title('Eccentricity')
end

sig={'Dr = 40%','Dr = 75%'};
figure(111)
legend(sig)

% Load

tests = [];
%tests = [1 2 4 5];

for t=1:numel(tests)
    test = tests(t);    
    S = Results{test};
    Device = S.Device;
    
    length_array = size(Device{1},2);
    crop_length = 2.5; % On both sides
    vp = [round(length_array*(crop_length/100)) ; length_array - round(length_array*(crop_length/100))];
    xv = linspace(crop_length,100-crop_length,vp(2)-vp(1)+1);       
    
    for s=1:numel(Device)
        Indexes = get_probe_geometry(Device{s});
        
        figure(30+t)
        subplot(1,3,1)
        hold on
        plot(xv,Indexes(1,vp(1):vp(2)))
        title('Cross section area')
        xlabel('Percentage of length') 
        ylabel('Area [mm2]') 
        
        subplot(1,3,2)
        hold on
        plot(xv,Indexes(2,vp(1):vp(2)))
        title('Major Axis')
        xlabel('Percentage of length') 
        ylabel('Axis length [mm]') 
        
        subplot(1,3,3)
        hold on
        plot(xv,Indexes(2,vp(1):vp(2))./Indexes(3,vp(1):vp(2)));
        title('Ratio between Axis')
        xlabel('Percentage of length')
            
    end
    
    legend('Step 1','Step 2','Step 3')
    sgtitle(['Effect of load - ' tests_list{test}])
end

% Single test L6
tests = [10 11 13];

for c=1:numel(tests)
    test = tests(c);    
    S = Results{test};
    Device = S.Device;
    
    length_array = size(Device{1},2);
    crop_length = 2.5; % On both sides
    vp = [round(length_array*(crop_length/100)) ; length_array - round(length_array*(crop_length/100))];
    xv = linspace(crop_length,100-crop_length,vp(2)-vp(1)+1);       
    
    s = 3;
         
    Indexes = get_probe_geometry(Device{s});
    
    figure(121)
    subplot(1,3,1)
    hold on
    plot(xv,Indexes(1,vp(1):vp(2)))
    xlim([0 100])
    %ylim([80 200])
    xlabel('Percentage of length') 
    ylabel('Area [mm2]') 
    title('Cross section area [mm2]')

    subplot(1,3,2)
    hold on
    plot(xv,Indexes(2,vp(1):vp(2)))
    xlim([0 100])
    %ylim([10 20])
    xlabel('Percentage of length') 
    ylabel('Length [mm]') 
    title('Major Axis length [mm]')
    
    subplot(1,3,3)
    hold on
    plot(xv,Indexes(2,vp(1):vp(2))./Indexes(3,vp(1):vp(2)));
    xlim([0 100])
    %ylim([1 1.2])
    xlabel('Percentage of length') 
    ylabel('Ratio of major/minor axis') 
    title('Eccentricity')
end

sig={'\sigma = 3.5 kPa','\sigma = 6 kPa','\sigma = 8.5 kPa'};
figure(121)
legend(sig)

%% Auxiliary Functions

function [Indexes] = get_probe_geometry(Im)
    %Ao = (pi/4)*10^2; % Area in mm2
    
    Im = permute(Im,[3 1 2]);
    nS = size(Im,3);
    
    Indexes = zeros(3,nS); % first row area second row shape
    
    for i=1:nS
        I = (Im(:,:,i))>0;
        stats = regionprops(I,'Area','MajorAxisLength','MinorAxisLength');
        if ~isempty(stats)           
            Indexes(:,i) = [ stats.Area*(0.07^2); stats.MajorAxisLength*0.07 ; stats.MinorAxisLength*0.07 ];
        else
            Indexes(:,i) = [NaN, NaN, NaN ];
        end
    end
end

function [table] = get_ledger(images_path,Disp_Folder,Strain_Folder)
%% 1. Get structure of the Tests
% find tests and their steps

srcFiles = dir([images_path '/*.tif']); % Weak check, to improve
table = cell(numel(srcFiles),8);
list = cell(numel(srcFiles),2);
for s = 1:numel(srcFiles)
   array = strsplit(srcFiles(s).name,'_');
   list(s,1:2)={array{2},array{3}(1:end-4)};
end

tests=unique(list(:,1)) % tests included in the results folder (disp)
line=0;

param = '_hws16_psr12_ns24_';

DispFiles = extractfield(dir([Disp_Folder '/' '*_hws16_psr12_ns24_*.tsv']),'name'); % Weak check, to improve
StrainFiles = extractfield(dir([Strain_Folder '/' '*_hws16_psr12_ns24_*.tif']),'name');

for t = 1 : numel(tests)  
    
    indexes = find(contains(list(:,1),tests{t}));

    for i = 1:numel(indexes)-1
        line = line+1;
        table(line,1:3)={tests{t},list{indexes(i),2},list{indexes(i+1),2}};
        
        % Disp 
        Bin1=contains(DispFiles,[table{line,1} '_' table{line,2} '_' table{line,3} param 'Bin1.tsv']);
        U4=contains(DispFiles,[table{line,1} '_' table{line,2} '_' table{line,3} param 'U4.tsv']);
        if any(U4)
            table{line,4}=2;
            table{line,6}=DispFiles{find(U4)};
        elseif any(Bin1)
            table{line,4}=1;
            table{line,6}=DispFiles{find(Bin1)};
        else
            table{line,4}=0;
        end
        
        rs=dlmread(fullfile(Disp_Folder,table{line,6}),'\t',1,18);
        table{line,8}=sum(and(rs(:,2)<0,rs(:,2)>-4))/sum(rs(:,2)>-4)*100;
        
        % Strain
        Bin1=contains(StrainFiles,[table{line,1} '_' table{line,2} '_' table{line,3} param 'Bin1' '-deviatoric-largeStrains.tif']);
        U4=contains(StrainFiles,[table{line,1} '_' table{line,2} '_' table{line,3} param 'U4' '-deviatoric-largeStrains.tif']);
        if any(U4)
            table{line,5}=2;  
            table{line,7}=StrainFiles{find(U4)};
        elseif any(Bin1)
            table{line,5}=1;
            table{line,7}=StrainFiles{find(Bin1)};
        else
            table{line,5}=0;
        end
        
        % Percentage of non convergence
        
    end
    
end

table(line+1:end,:)=[];

end

function [I]=load_tif_stack(path)
    InfoImage=imfinfo(path);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
    I=zeros(nImage,mImage,NumberImages);

    TifLink = Tiff(path, 'r');
    for i=1:NumberImages
       TifLink.setDirectory(i);
       I(:,:,i)=TifLink.read();
    end
    TifLink.close();
end

function [Top,Front,Side] = get_slices(I)
   
    slice=round(size(I,1)/2);
    
    Top=I(:,:,round(size(I,3)/2));
    Side=flipud(permute(I(slice,:,:),[3 2 1])); % Side   
    Front=flipud(permute(I(:,slice,:),[3 1 2])); % Front

end

function [Disp, Strain]=get_Images(disp_path,strain_path,DevicePath)
    
    %Device Mask
    SoilStruct=load(fullfile(DevicePath,'Soil.mat'));
    SDM=SoilStruct.SDM;
    Mask=SoilStruct.SCM;
    Device=zeros([size(Mask) 1200]);
    Device(SDM)=1;
    
    % Displacement
    if contains(disp_path,'Bin')
        nSlices=1200;  
    else
        nSlices=1000;
        Device(:,:,1001:end)=[];
        Device=flipud(Device);
    end
    
    suffix_Stat='-SubPixReturnStat.tif';
    Stat=load_tif_stack([disp_path suffix_Stat]);  
    
    if size(Stat,3)>100
        error
    end
    
    Device = imresize3(Device,size(Stat));   
  
    suffix_disp={'-Zdisp.tif','-Ydisp.tif','-Xdisp.tif'};
    Disp=zeros([size(Stat) 4]);
    Im = zeros(size(Stat));
    
    for i=1:3
        I=load_tif_stack([disp_path suffix_disp{i}]);
        I(Stat<1.5)=Inf; % Non-convergence
        %I(Device>0.05)=Inf; % Mask
        Disp(:,:,:,i)=I;
        Im = I.^2 + Im;
    end

    Im=Im.^(1/2);
    Im(Stat<1.5)=Inf; % Non-convergence
    %Im(Device>0.05)=Inf; % Mask
    Disp(:,:,:,4)=Im;
       
    % Strain   
    
    suffix_strain={'deviatoric','volumetric','xx','xy','xz','yx','yy','yz','zx','zy','zz','-'};
    Strain=zeros([(size(Stat)-1) (numel(suffix_strain)-1)]);

    Stat = imresize3(Stat,size(Stat)-1);
    Device = imresize3(Device,size(Device)-1);
    
    for i=1:numel(suffix_strain)-1
        
        I=load_tif_stack(strain_path);
        I(Stat<1.5)=Inf; % Non-convergence
        %I(Device>0.05)=Inf; % Mask
        Strain(:,:,:,i)=I;

        strain_path = strrep(strain_path,suffix_strain{i},suffix_strain{i+1});
    end
    
end

function [names]=GetSubfolders(path)
        d = dir(path);
        isub = [d(:).isdir]; %# returns logical vector
        names = {d(isub).name}';
        names(ismember(names,{'.','..'})) = [];
end

function []=SaveMultiTiff(path,array)
    imwrite(array(:,:,1), path)
    for i=2:length(array(1,1,:))
        imwrite(array(:,:,i), path, 'writemode', 'append')
    end                       
end

function [MP] = GetPlainStrain(Strain)
% MP: [ Norm, Loss, eigenvectors (9), eigenvalues(3) max_shear]
% Strain: {'deviatoric','volumetric','xx','xy','xz','yx','yy','yz','zx','zy','zz'}
% Gets the ratio between the euclidean norms of the 2D and 3D tensors at
% every point of every test

MP = cell(size(Strain));

for t=1:size(MP,1) % Loop Tests
    for s=1:size(MP,2) % Loop Steps
        
        S = Strain{t,s};
        
        if isempty(S); continue; end
        
        Props = nan([size(S,1) size(S,2) size(S,3) 15]);
        
        for i=1:size(S,1)
            for j=1:size(S,2)
                for k=1:size(S,3)
                    v = permute(S(i,j,k,:),[4,3,2,1]);
                    % [ xx xz xy ; zx zz zy ; yx yz yy ]
                    %tens = [v(3) v(5) v(4); v(9) v(11) v(10); v(6) v(8) v(7)]; 
                    tens = [v(7) v(8) v(6); v(10) v(11) v(9); v(4) v(5) v(3)]; 
                    
                    if all(isfinite(tens(:)))
                        [V,D] = eig(tens); [d,ind] = sort(diag(D),'descend'); Vs = V(:,ind);
                        Props(i,j,k,1)= norm(tens,2); % Norm Tensor
                        Props(i,j,k,2)= 100*(1-norm(tens(1:2,1:2),2)/norm(tens,2)); % Loss                 
                        Props(i,j,k,3:end-1) = [Vs(:)' d'];
                        Props(i,j,k,end) = abs(d(1)-d(3));
                        
                    end
                end
            end
        end
        
        MP{t,s} = Props;
    end
end

end

function [m] = getmaxmin(X)

    A=X(:,:,:,1); % Deviatoric (positive)
    maxD = max(A(isfinite(A)),[],'all');

    A=X(:,:,:,2); % Volumetric (positive and negative)
    maxV = max(A(isfinite(A)),[],'all');
    minV = min(A(isfinite(A)),[],'all');

    m=[maxD minV maxV];

end


function []=PlotPrincStrainOrientation(Strain,Tests,step)

    try 
        load('streamlines.mat')
    catch
        streamlines = cell(13,3,2);
    end
    
    %divMap = [brewermap(101,'RdBu')]; % Divergent Colormap (- to +)
    seqMap = brewermap(101,'YlGn'); % Sequential Colormap (0 to inf)

    components = [ 1 2 ; 3 2 ; 3 1];
    sizes = [100 70 100];
    colors = {[64 64 64]./255,[153 0 0]./255};
       
    ntests = numel(Tests);
    FIG = figure;
    
    for ti = 1:ntests
        
        test = Tests(ti);       
        A = Strain{test,step}(:,:,:,12:end);

        ix = 100/size(A,1);
        xc = ix/2:ix:100-ix/2;
        yc = xc;
        iy = 70/40;%size(A,3);
        zc = 1.*(iy/2:iy:70-ix/2);

        [Grids{1,1}, Grids{1,2}] = meshgrid(xc,zc); % Front
        [Grids{2,1}, Grids{2,2}] = meshgrid(yc,zc); % Side
        [Grids{3,1}, Grids{3,2}] = meshgrid(xc,yc); % Top

        iivs = [1 3];
        for  iiv=1:2 % Major and minor principal strain
            iv = iivs(iiv);
            %P = cat(4,A(:,:,:,iv*3:iv*3+2).*A(:,:,:,11+iv),abs(A(:,:,:,12)-A(:,:,:,14))); % Principal Strain Array
            P = cat(4,A(:,:,:,iv*3:iv*3+2),abs(A(:,:,:,12)-A(:,:,:,14))); % Principal Strain Array      

            Slices{1} = flipud(permute(P(:,round(size(P,2)/2),1:40,:),[3 1 4 2])); % Front
            Slices{2} = flipud(permute(P(round(size(P,2)/2),:,1:40,:),[3 2 4 1])); % Side  
            Slices{3} = P(:,:,round(size(P,3)/2),:); % Top

            for view=1:3 % Front, Side and Top
                u = Slices{view}(:,:,components(view,1)); v = Slices{view}(:,:,components(view,2));
                shear = imresize(Slices{view}(:,:,4),[sizes(components(view,2)) sizes(components(view,1))]);

                Xc = Grids{view,1}; Yc = Grids{view,2};
                ys = sign(Yc-50.1); xs = sign(Xc-50.1); 

                if iv==1
                    sy = 1;
%                     [sxv,syv] = meshgrid([-2 2]+50,1:2:70);
%                     [sxh,syh] = meshgrid(0:2:100,[-2 5]+50);
                else
                    sy = -1;           
%                     [sxv,syv] = meshgrid([2 98],1:2:70);
%                     [sxh,syh] = meshgrid(0:2:100,[1 69]);
                end

                b = (sign(u)==xs)+2.*(sy.*sign(v)==ys);
                uc = u ; vc = v ;
                %uc(b==3) = uc(b==3); vc(b==3) = vc(b==3); % Good -> No change
                uc(b==0) = -1.*uc(b==0); vc(b==0) = -1.*vc(b==0); % Change both (opposite)
                %uc(b==2) = -1.*uc(b==2); vc(b==2) = -1.*vc(b==2); % Change both (opposite)
                uc(or(b==1,b==2)) = nan; vc(or(b==1,b==2)) = nan; % Change trouble-some (nan)
                %uc(b==2) = nan; vc(b==2) = nan; % Change trouble-some (nan)

                ucc = inpaint_nans(uc); vcc = inpaint_nans(vc); % Interpolate vector components

                if ~isempty(streamlines{test,view,iiv})
                    q = figure;
                    quiver(Xc,Yc,ucc,vcc)
                    axis image; set(gca, 'YDir','reverse')
                    streamlines{test,view,iiv} = readPoints(q,Xc,Yc,ucc,vcc);
                    close(q);
                end
                
                save('streamlines.mat','streamlines');
                dhfjg=1;

                %contour = drawpolyline;
                %x=contour.Position(:,1); y=90-contour.Position(:,2);
                
                figure(FIG)
                f = subplot(ntests,3,(ti-1)*3+view); 
                    hold on
                    h = streamline(stream2(Xc,Yc,ucc,vcc,streamlines{test,view,iiv}(:,1),streamlines{test,view,iiv}(:,2)));
                    set(h,'Color',colors{iiv})
                    axis image; set(gca,'YDir','reverse')

    %             f = subplot(1,3,view);   
    %             %imagesc(shear);
    %             %colormap(f,seqMap);
    %             axis image ; hold on;
    %             q=quiver(Grids{view,1},Grids{view,2},u,v,1) ;
    %             q.Color = colors{iv}; %q.ShowArrowHead = 'off';
    %             set(gca, 'YDir','reverse')
    %             dhd=1;

            end

        end

    end
    
    save('streamlines.mat','streamlines');
end


function   [xc,yc,R,a] = circfit(x,y)
%
%   [xc yx R] = circfit(x,y)
%
%   fits a circle  in x,y plane in a more accurate
%   (less prone to ill condition )
%  procedure than circfit2 but using more memory
%  x,y are column vector where (x(i),y(i)) is a measured point
%
%  result is center point (yc,xc) and radius R
%  an optional output is the vector of coeficient a
% describing the circle's equation
%
%   x^2+y^2+a(1)*x+a(2)*y+a(3)=0
%
%  By:  Izhak bucher 25/oct /1991, 
    x=x(:); y=y(:);
   a=[x y ones(size(x))]\[-(x.^2+y.^2)];
   xc = -.5*a(1);
   yc = -.5*a(2);
   R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
   
end

function [a,b,c,rsq] = fit_catenary(x,y)

    Eqn = 'a*cosh((x-b)/a)+c';
    [f1,gn] = fit(x,y,Eqn,'StartPoint',[25 50 10]);
    coeff = coeffvalues(f1);
    a = coeff(1); b = coeff(2); c = coeff(3);
    rsq = gn.rsquare;
    
end



function pts = readPoints(fig,Xc,Yc,ucc,vcc)
%readPoints   Read manually-defined points on current image, streamlines
%   appear as points are clicked on the quiver plot until the right button
%   is clicked by the user.
% 

    pts = zeros(2, 0);
    k = 0;
    figure(fig)
    hold on;           % and keep it there while we plot

    while 1
        [xi, yi, but] = ginput(1);      % get a point

        if ~isequal(but, 1)             % stop if not button 1
            break
        end

        k = k + 1;
        pts(1,k) = xi; pts(2,k) = yi;

        streamline(stream2(Xc,Yc,ucc,vcc,xi,yi));

    end

    hold off;
    pts = pts';
end
