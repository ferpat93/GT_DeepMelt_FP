
% Check what I got
close all;
Folders={'2ND-WEEK-'};
root='G:\';
images_path = fullfile(root,'DIC','source-images','Under4');
Disp_Folder=fullfile(root,'DIC','results','displacement');
Strain_Folder=fullfile(root,'DIC','results','strain');
segmentation_folder = 'G:\2ND-WEEK-OP2';
Scan_Folder = fullfile(root,strcat(Folders{1},'SC'),filesep);    
    
Thresholds = xlsread('Porosity_thresholds.xlsx',1,'B1:D13');

hws = 12; % 3/4 the value used in DIC.
ns = 1.5*hws; % Innacurate volume measurement if not equal to 2

%{ 
Types of value:
 0: Not computed yet
 1: Computed using big images "Bin1"
 2: Computed using reduced slices 'Under4'
%}

table = get_ledger(images_path,Disp_Folder,Strain_Folder); % Ledger of DIC results
tests_list = unique(table(:,1));

%% Get Results

Disp = cell(numel(tests_list),3);
Strain = cell(numel(tests_list),3);
%PV
Porosity = cell(numel(tests_list),4);
Volume = zeros(numel(tests_list),4);
Device = cell(numel(tests_list),4);

for t = 1:numel(tests_list)
    
    disp(tests_list{t})
    
    % Get Displacement and Strain
    indexes = find(contains(table(:,1),tests_list{t}));
    for i=1:numel(indexes)   
        if and(table{indexes(i),4}>0,table{indexes(i),5}>0)
            DeviceFolder = fullfile(segmentation_folder,table{indexes(i),1},table{indexes(i),2});
            disp_path = fullfile(Disp_Folder,table{indexes(i),6}(1:end-4));
            strain_path = fullfile(Strain_Folder,table{indexes(i),7});
            [Disp{t,i}, Strain{t,i}]=get_Images(disp_path,strain_path,DeviceFolder);
        end
    end 
    
    % Get Porosity      
    TT=Thresholds(t,:);
    nameSteps=GetSubfolders(fullfile(Scan_Folder,tests_list{t}));

    for s=1:numel(nameSteps) % Iterates Though Steps of a Test
        disp(nameSteps{s}); 
        ImagesFolder=fullfile(Scan_Folder,tests_list{t},nameSteps{s},'SlicesY');
        ResultsFolder=fullfile(segmentation_folder,tests_list{t},nameSteps{s});

        if s==1
            D=load(fullfile(ResultsFolder,'Device.mat'));
            DB = uint16(D.DB);
        end
        Soil=load(fullfile(ResultsFolder,'Soil.mat'));

        [Stack,XCidx,Centroid]=GetVolume_OtherDirection(ImagesFolder,Soil,TT(3),DB);
        Stack(Soil.SDM)=0;
        Device{t,s}=XCidx;
        [Porosity{t,s},Volume(t,s)] = GetPorosityMap(Stack,TT(1:2),ns,hws);

    end

end

save('Data.mat','Porosity','Disp','Strain','Volume','Device');

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

tests=unique(list(:,1)); % tests included in the results folder (disp)
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

function [PS,MS] = GetPlainStrain(Strain)
% {'deviatoric','volumetric','xx','xy','xz','yx','yy','yz','zx','zy','zz'}
% Gets the ratio between the euclidean norms of the 2D and 3D tensors at
% every point of every test

PS = cell(size(Strain));
MS = cell(size(Strain));

for t=1:size(PS,1) % Loop Tests
    for s=1:size(PS,2) % Loop Steps
        
        S = Strain{t,s};
        
        if isempty(S); continue; end
        
        loss = zeros([size(S,1) size(S,2) size(S,3)]);
        mag = zeros([size(S,1) size(S,2) size(S,3)]);
        
        for i=1:size(S,1)
            for j=1:size(S,2)
                for k=1:size(S,3)
                    v = permute(S(i,j,k,:),[4,3,2,1]);
                    % [ xx xz xy ; zx zz zy ; yx yz yy ]
                    %tens = [v(3) v(5) v(4); v(9) v(11) v(10); v(6) v(8) v(7)]; 
                    tens = [v(7) v(8) v(6); v(10) v(11) v(9); v(4) v(5) v(3)]; 
                    loss(i,j,k)= 100*(1-norm(tens(1:2,1:2),2)/norm(tens,2));
                    mag(i,j,k)= norm(tens,2);
                    
                end
            end
        end
        
        loss(isnan(loss)) = nan;
        mag(isinf(mag)) = nan;
        PS{t,s} = loss;
        MS{t,s} = mag;
    end
end

end
