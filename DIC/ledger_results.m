
% Check what I got
close all;
root = 'G:\DIC\';
images_path = fullfile(root,'source-images','Under4');

Disp_Folder=fullfile(root,'results','displacement');
Strain_Folder=fullfile(root,'results','strain');
segmentation_folder = 'G:\2ND-WEEK-OP2';

%{ 
Types of value:
 0: Not computed yet
 1: Computed using big images "Bin1"
 2: Computed using reduced slices 'Under4'
%}

table = get_ledger(images_path,Disp_Folder,Strain_Folder);
tests_list = unique(table(:,1));
%% Get Results

GetResults = 1;
if ~GetResults
    load('Results.mat');
else
    Results = cell(numel(tests_list),1);
    for t = 1:numel(tests_list)
        tests_list{t}
        indexes = find(contains(table(:,1),tests_list{t}));
        S = struct;
        S.numSteps = numel(indexes)+1;
        S.Indexes = indexes;

        Disp = cell(numel(indexes),1);
        Strain = cell(numel(indexes),1);

        for i=1:numel(indexes)   
            if and(table{indexes(i),4}>0,table{indexes(i),5}>0)
                DeviceFolder = fullfile(segmentation_folder,table{indexes(i),1},table{indexes(i),2});
                disp_path = fullfile(Disp_Folder,table{indexes(i),6}(1:end-4));
                strain_path = fullfile(Strain_Folder,table{indexes(i),7});
                [Disp{i}, Strain{i}]=get_Images(disp_path,strain_path,DeviceFolder);
            end
        end 

        Device = cell(numel(indexes)+1,1);
        steps = {table{indexes,2}, table{indexes(end),3}};
        TestCell = load(fullfile(segmentation_folder,strcat(tests_list{t},'.mat')));
        TestCell = TestCell.TestCell;

        for i=1:numel(steps)   
            DevicePath = fullfile(segmentation_folder,tests_list{t},steps{i},'SmoothDevice.tif');
            Device{i} = load_tif_stack(DevicePath);
            %Device{i} = TestCell{i}.CenteredCloud;
        end 

        S.Disp = Disp;
        S.Strain = Strain;
        S.Device = Device;
        Results{t} = S;
    end
end


%% Plain Strain Soil
 tests = [4 3 8];
% tests = [];

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