%% Device Analysis

close all

%% Lines to use tight subplots
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.065 0.065], [0.2 0.05], [0.075 0.04]); % subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
% if ~make_it_tight,  clear subplot;  end

load('Data.mat');
%load('StrainData.mat');

TestVariables = [ [2 4 6 2 4 6 2 4 2 6 6 4 6]' [1 1 1 1 1 1 0 0 0 0 0 0 0]' [0 0 0 1 1 1 0 0 1 0 1 1 2]']; % Length (2,4,6), density (0,1), load(0,1,2)

%% Device Shape

% 3 plots - one per length
% Area - AR - Centroid

% [yC, Centroid(2), major-minor axis, orientation, area, solidity]
steps = [2 3];
InitialArea=(3.1416/4)*(10.5)^2;

for l=8:2:6
    
    idx = find(TestVariables(:,1)==l);
    
    for i=1:numel(idx)
        
        if isempty(Device{idx(i),1}); continue; end
        
        Do = Device{idx(i),1};
        Co = mean(Do(:,3));
        bounds = [min(Do(Do(:,1)>0,1)) max(Do(:,1))];
        center = mean(bounds);
        
        figure;
        for s=1:numel(steps)
            si=steps(s);
            K = 10;
            
            D = Device{idx(i),si}; D(D(:,8)==0,:)=[];
            Shape = movmean(D(:,5)./D(:,4),K);
            Area = movmean(D(:,7)/InitialArea,K);
            Centroid = movmean(D(:,3)-Co,K);
            xv = (D(:,1)-center)*0.07;
            
            %subplot(1,3,round(l/2)); hold on % Centroid
            subplot(1,3,1); 
            hold on % Centroid
            plot(xv,Centroid)
            %ylim([-5 5]) ; 
            xlim([-5*l 5*l])
            
            subplot(1,3,2); 
            hold on % Area
            plot(xv,Area)
            ylim([1 2.6]); 
            xlim([-5*l 5*l])
            
            subplot(1,3,3);
            hold on % Shape
            plot(xv,Shape)
            %ylim([0.7 1]); 
            xlim([-5*l 5*l])
            
        end
        
    end
    
end


%% Old Code - Shape and Extent

load_indexes = true;

if load_indexes
    load('indexes_device.mat');
else

    load('Old_Device.mat');

    %tests = [];
    tests = 1:13;
    indexes_device = cell(13,1);

    % Per test
    for t=1:numel(tests)
        test = tests(t);    
        S = Results{test};
        Device = S.Device;

        length_array = size(Device{1},2);
        crop_length = 2.5; % On both sides
        vp = [round(length_array*(crop_length/100)) ; length_array - round(length_array*(crop_length/100))];
        xv = linspace(crop_length,100-crop_length,vp(2)-vp(1)+1);       

        s = 3; % Only for second step

        %for s=2:numel(Device)
            DeviceStep = Device{s};
            Indexes = get_old_probe_geometry(Device{s});

            indexes_device{test} = [xv; Indexes(:,vp(1):vp(2))];

    end
end

tests_to_plot = [ 1 2 3 ; 10 11 nan ; 9 4 nan];
legends{3} = {'Dr = 40%','Dr = 75%'};
legends{2} = {'\sigma = 3.5 kPa','\sigma = 6 kPa'}; % ,'\sigma = 8.5 kPa'
%legends{1} = {'L = 2 cm','L = 4 cm','L = 6 cm'};
legends{1} = {'L = 2D','L = 4D','L = 6D'};

titles = {'Device Length','Vertical stress','Soil Density'};

colorlines = {'r','k','b'};
linetyp = {'-.','-','--'};



figure;
for pi = 1:3
       
    for si = 1:sum(~isnan(tests_to_plot(pi,:)))
        
        ddd = indexes_device{tests_to_plot(pi,si)};
           
        subplot(2,3,pi); hold on
        plot(ddd(1,:),ddd(2,:),[colorlines{si} linetyp{si}])
        ylim([1.25 2.5]); grid on
        yticks(1.25:0.25:2.5)
        xtickformat('%g%%');
        xticks(0:20:100)
        title(titles{pi})
        xlabel('Device length')
        ylabel('Cross-section area')
        
        subplot(2,3,pi+3); hold on
        plot(ddd(1,:),ddd(4,:),[colorlines{si} linetyp{si}])  
        ylim([0 20]); grid on
        ytickformat('%g%%');
        xtickformat('%g%%');
        xticks(0:20:100)
        yticks(0:5:20)
        xlabel('Device length')
        ylabel('Eccentricity')
        
    end
    
    legend(legends{pi},'Location','south')

end

tableI = zeros(numel(indexes_device),4);

for i=1:numel(indexes_device)
    
    di = indexes_device{i}';
    
    ma = mean(di(:,2));
    ms = mean(di(:,4));
    
    ra = (max(di(20:end-20,2))/ma)-1;
    rs = (max(di(20:end-20,4))/ms)-1;
    
    tableI(i,:) = [ma ra ms rs];
end


%% AUXILIAR FUNCTIONS

% 

function [Indexes] = get_old_probe_geometry(Im)
    %Ao = (pi/4)*10^2; % Area in mm2
    
    Im = permute(Im,[3 1 2]);
    nS = size(Im,3);
    InitialArea=(3.1416/4)*(10.5)^2;
    
    Indexes = zeros(3,nS); % first row area second row shape
    
    for i=1:nS
        I = (Im(:,:,i))>0;
        stats = regionprops(I,'Area','MajorAxisLength','MinorAxisLength');
        if ~isempty(stats)           
            Indexes(:,i) = [ stats.Area.*(0.07^2)./InitialArea; stats.MajorAxisLength*0.07 ; 100.*((stats.MajorAxisLength/stats.MinorAxisLength)-1) ];
        else
            Indexes(:,i) = [NaN, NaN, NaN ];
        end
    end
end