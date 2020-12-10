
close all

%% Lines to use tight subplots
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.065 0.065], [0.2 0.05], [0.075 0.04]); % subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
if ~make_it_tight,  clear subplot;  end

load('Data.mat');
load('StrainData.mat');

TestVariables = [ [2 4 6 2 4 6 2 4 2 6 6 4 6]' [1 1 1 1 1 1 0 0 0 0 0 0 0]' [0 0 0 1 1 1 0 0 1 0 1 1 2]']; % Length (2,4,6), density (0,1), load(0,1,2)

%% Porosity

divMap = [brewermap(101,'RdBu')];
seqMap = [brewermap(101,'YlGn')];

for t=8:size(Porosity,1)
    
    figure(t)
    
    P{1} = Porosity{t,1};
    P1 = Porosity{t,2};
    P2 = Porosity{t,3};
    
    P{2} = P1-P{1};
    P{3} = P2-P{1};
    
    for i=1:3
           
        [~,Front,Side] = get_slices(P{i});  % Volumetric
        %maglim = max(abs(extremes(t,2,step)),extremes(t,3,step)); lim = [-maglim maglim];
        maglim = max(abs(prctile([Front(:);Side(:)],5)),prctile([Front(:);Side(:)],95));
        if i==1; lim = [0 100]; CM =  seqMap; else; lim = [-maglim maglim]; CM =  divMap;  end
        ax(i)=subplot(2,3,i); imagesc(Front,lim); colormap(ax(i),CM);axis image ; axis off;
        ax(i+3)=subplot(2,3,i+3); imagesc(Side,lim); colormap(ax(i+3),CM);axis image ; axis off;
        colorbar('southoutside')

    end
end 



%% Volumetric vs Deviatoric Strain
%%% -> Compressive volumetric strains are negative
% NOT USED FROM NOW BC EXTREMES ARE NOT IN SECTIONS
%A = cellfun(@(x) getmaxmin(x),Strain(:,1:2),'UniformOutput',false);
%extremes(:,:,1)=cell2mat(A(:,1));
%extremes(:,:,2)=cell2mat(A(:,2)); % Max Deviatoric - Min Volumetric - Max Volumetric

% Plot to check strain maps

step=2;
divMap=[[0 0 0]; brewermap(101,'RdBu')];
seqMap=[brewermap(101,'YlGn');[0 0 0]];

for t=8:size(Strain,1)
    figure(t+50)
    
    %Deviatoric
    [~,Front_D,Side_D] = get_slices(Strain{t,step}(:,:,:,1));  % Deviatoric
    %lim = [0 extremes(t,1,step)*1.01];
    lim = [ 0 1.01*prctile([Front_D(isfinite(Front_D));Side_D(isfinite(Side_D))],97.5)];
    ax(1)=subplot(2,3,1); imagesc(Front_D,lim); colormap(ax(1),seqMap);axis image ; axis off;
    ax(4)=subplot(2,3,4); imagesc(Side_D,lim); colormap(ax(4),seqMap);axis image ; axis off;
    colorbar('southoutside')
    
    % Volumetric
    [~,Front_V,Side_V] = get_slices(Strain{t,step}(:,:,:,2));  % Volumetric
    %maglim = max(abs(extremes(t,2,step)),extremes(t,3,step)); lim = [-maglim maglim];
    maglim = max(abs(prctile([Front_V(isfinite(Front_V));Side_V(isfinite(Side_V))],5)),prctile([Front_V(isfinite(Front_V));Side_V(isfinite(Side_V))],95));
    lim = [-maglim maglim];
    ax(2)=subplot(2,3,2); imagesc(-1*Front_V,lim); colormap(ax(2),divMap);axis image ; axis off;
    ax(5)=subplot(2,3,5); imagesc(-1*Side_V,lim); colormap(ax(5),divMap);axis image ; axis off;
    colorbar('southoutside')
    
    % Ratio
    FR = -1*Front_D./Front_V;
    SR = -1*Side_D./Side_V;
    
    maglim = max(abs(prctile(FR(:),10)), abs(prctile(FR(:),90))); lim = [-maglim maglim];
    ax(3)=subplot(2,3,3); imagesc(FR,lim); colormap(ax(3),divMap);axis image ; axis off;
    maglim = max(abs(prctile(SR(:),10)), abs(prctile(SR(:),90))); lim = [-maglim maglim];
    ax(6)=subplot(2,3,6); imagesc(SR,lim); colormap(ax(6),divMap);axis image ; axis off;
    colorbar('southoutside')
        
    %[~,Front_M,Side_M] = get_slices(MP{t,step}(:,:,:,1));  % Magnitude
    figure;
    imagesc(imresize(Front_D,round(1.75*size(Front_D))));
    h = drawfreehand;
    cc = createMask(h);
    close;
end 

%% Plain Strain Soil
% (Probe axis is along Y - 2nd dimension)
% Norm difference (2D vs. 3D sensor), figures alongside percentages of
% variation.

v=0:(100/58):100;
pPS = zeros(13,2);
testsToPlot = [1 2 3];
ImagesToPlot = cell(2,numel(testsToPlot),2);
maxStrain = [0 0]; steptoplot = 1;

colorlines = {'r','k','b'};
linetyp = {'-.','-'};

SM = cell(size(MP)); PS = cell(size(MP));
for ii=1:numel(MP)
    if ~isempty(MP{ii})
        SM{ii} = MP{ii}(:,:,:,1);
        PS{ii} = MP{ii}(:,:,:,2);
    end
end

for l=2:2:6
    
    idx = find(TestVariables(:,1)==l);
    intv = (-l*5)-5:1.25:(l*5)+5;
    
    for i=1:numel(idx)
        
        D = Device{idx(i),1}; center = 0.07*(D(1,3)+D(1,2))/2;
        vc = v - center;
            
        for s=1:2
            
            meanLoss = nansum(PS{idx(i),s}.*SM{idx(i),s},[2,3])./nansum(SM{idx(i),s},[2,3]);
            intLoss = interp1(vc,meanLoss,intv); 
            
            figure(s)
            subplot(1,3,l/2)
            hold on
            plot(intv,intLoss,[colorlines{TestVariables(idx(i),3)+1} linetyp{TestVariables(idx(i),2)+1}])
            
            pPS(idx(i),s) = 100*(sum(abs(intLoss)<=5)/sum(abs(intv)<=l*5));
            
        end
        
        ip = find(testsToPlot==idx(i));
        if ~isempty(ip)   
            mag = SM{idx(i),steptoplot}; psi = PS{idx(i),steptoplot}; 
            maxStrain(1) = max(maxStrain(1),max(mag(~isinf(mag)))); maxStrain(2) = max(maxStrain(2),max(psi(~isinf(psi))));
            mag(isnan(mag)) = inf;   psi(isnan(psi)) = inf;       
            [~,ImagesToPlot{1,ip,1},ImagesToPlot{2,ip,1}] = get_slices(mag);
            [~,ImagesToPlot{1,ip,2},ImagesToPlot{2,ip,2}] = get_slices(psi); 
            
        end

    end
    
    figure(1); subplot(1,3,l/2); xlim([min(intv) max(intv)]); ylim([0 60]); xlabel('Distance from axis center [mm]'); ylabel('Average Plane Strain Deviation [%]'); title(strcat('L =',{' '},num2str(l),'D')); grid on
    figure(2); subplot(1,3,l/2); xlim([min(intv) max(intv)]); ylim([0 60]); xlabel('Distance from axis center [mm]'); ylabel('Average Plane Strain Deviation [%]'); title(strcat('L =',{' '},num2str(l),'D')); grid on
end

for s=1:2; figure(s); subplot(1,3,3); legend('Dense - s1','Dense - s2','Loose - s1','Loose - s2','Loose - s3'); end

% Plot Strain Magnitude and Loss

divMap=[brewermap(100,'Oranges');[0 0 0]];
maxStrain = [0.1 100];
for f = 1:2
    for li=1:3
        figure(20+f)
        subplot(2,3,li); imagesc(ImagesToPlot{1,li,f}(1:40,:),[0 maxStrain(f)]);  axis image ; axis off;
        title(strcat('L = ',num2str(2*li),'D'))
        subplot(2,3,li+3); imagesc(ImagesToPlot{2,li,f}(1:40,:),[0 maxStrain(f)]); axis image ; axis off;
        colormap(divMap)       
    end 
    %colorbar('eastoutside')
    %h(1)=text(0.05,0.7,'Transversal Section'); h(2)=text(-2.65,0.2,'Longitudinal Section');
    %h(1) = annotation('textbox',[.05 .7 .1 .2],'String','Text outside the axes','EdgeColor','none','FitBoxToText','on');
    %set(h(1),'Rotation',90); set(h(2),'Rotation',90);
end

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
% [ Norm, Loss, eigenvectors (9), eigenvalues(3)]
% {'deviatoric','volumetric','xx','xy','xz','yx','yy','yz','zx','zy','zz'}
% Gets the ratio between the euclidean norms of the 2D and 3D tensors at
% every point of every test

MP = cell(size(Strain));

for t=1:size(MP,1) % Loop Tests
    for s=1:size(MP,2) % Loop Steps
        
        S = Strain{t,s};
        
        if isempty(S); continue; end
        
        Props = nan([size(S,1) size(S,2) size(S,3) 14]);
        
        for i=1:size(S,1)
            for j=1:size(S,2)
                for k=1:size(S,3)
                    v = permute(S(i,j,k,:),[4,3,2,1]);
                    % [ xx xz xy ; zx zz zy ; yx yz yy ]
                    %tens = [v(3) v(5) v(4); v(9) v(11) v(10); v(6) v(8) v(7)]; 
                    tens = [v(7) v(8) v(6); v(10) v(11) v(9); v(4) v(5) v(3)]; 
                    
                    if all(isfinite(tens(:)))
                        [V,D] = eig(tens); [d,ind] = sort(diag(D)); Vs = V(:,ind);
                        Props(i,j,k,1)= norm(tens,2); % Norm Tensor
                        Props(i,j,k,2)= 100*(1-norm(tens(1:2,1:2),2)/norm(tens,2)); % Loss                 
                        Props(i,j,k,3:end) = [Vs(:)' d'];
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

