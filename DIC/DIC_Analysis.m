% Color Coded Displacement analysis

root='G:';
Disp_Folder=fullfile(root,'DIC','results','displacement');
Strain_Folder=fullfile(root,'DIC','results','strain');

%% 1. Get structure of the folder
% find tests within folder, needed to loop through them

srcFiles = dir([Disp_Folder '/' '*ns24_Bin1.tsv']); % Weak check, to improve
files = cellfun(@(x) x(1:end-4),extractfield(srcFiles,'name'),'UniformOutput',false)'; % list of filenames
%split_files = cellfun(@(x) strsplit(x,'_'),files,'UniformOutput',false);
files_comp={};
for s = 1:numel(srcFiles)
   files_comp=[files_comp; strsplit(srcFiles(s).name,'_')]; % file components: test, step, hws,...
end
tests=unique(files_comp(:,2)); % tests included in the results folder (disp)

%% 2. Loop through tests,

% Components:
comp={'Zdisp','Xdisp','Ydisp','SubPixReturnStat'};
Images=cell(numel(tests),1);

% Loop over tests
for f=1:3  %f=1:numel(tests)
    index_tests=find(strcmp(files_comp(:,2),tests{f})); %find completed steps of test    
    seq_images=cell(numel(index_tests),1);
    
    for t=1:numel(index_tests)
        seq_images{t}=get_richImage(Disp_Folder,files{index_tests(t)});
    end  
    Images{f}=seq_images;
end


%% 3. Color Mapping and plotting

% Images -> Cell array of cell arrays: as many internal cells as tests,
% each test cell has n steps

%{ 

% Displacement (3 Components)
RGB COLOR CODING
    Red: Z disp
    Green: X disp
    Blue: Y disp
    
** White: not convergence (or max disp) **

%}

res=70/1000; % resolution of scan in mm

% Loop over tests
for f=1:3  %f=1:numel(tests) % loop over tests
    maxTest=max(cellfun(@(x) max(abs(x(:))),Images{f}));
    CK=255/maxTest; % Constant to map displacement to color
    
    for t=1:numel(Images{f}) % loop over steps
        [XZ, YZ, XY] = map_to_RGB(Images{f}{t},CK);
        [XZm, YZm, XYm] = magnitude_map(Images{f}{t});
        
        HI = histo_Axis(Images{f}{t});        
        S = {XZ,YZ,XY,XZm,YZm,XYm};
        plot_group(S,HI,[tests{f} ' - Step ',num2str(t)])    
      
    end  

end





%% Auxiliar Functions

function [I]=get_richImage(folder,name_path)

    % Components:
    comp={'Zdisp','Xdisp','Ydisp','SubPixReturnStat'};
    RS=load_tif_stack(fullfile(folder,[name_path '-' comp{numel(comp)} '.tif']));
    I=zeros(size(RS,1),size(RS,2),size(RS,3),4);
    I(:,:,:,4)=RS; 
    
    for i=1:numel(comp)-1 % Loads each disp component and masks non-convergent voxels
        Il=load_tif_stack(fullfile(folder,[name_path '-' comp{i} '.tif']));
        %Il(RS<0)=-Inf;
        I(:,:,:,i)=Il;
        % Disp(Stat<-3)=Inf; % For masked region
    end

    
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

function [XZ, YZ, XY]=map_to_RGB(I,C)
    
    % I is a 4D image {Z, X, Y, RS}
    
    %{ 
    % Displacement (3 Components)
    RGB COLOR CODING
        Red: Z disp
        Green: X disp
        Blue: Y disp
    ** White: not convergence (or max disp) **
    %}

    SI = zeros(size(I,1),size(I,2),size(I,3),3,'uint8'); % RGB matrix
    
    for c=1:3 % Loop over Components
        scaled_I = round(abs(I(:,:,:,c)).*C); % scale with given constant
        scaled_I(I(:,:,:,4)<0)= 255; % set non-conv to white (max value)
        SI(:,:,:,c)=scaled_I; % store    
    end
    
    scale=10;
    % extract slices
    YZ = imresize(flipud(permute(SI(round(size(SI,1)/2),:,:,:),[3 2 4 1])),scale); % Side
    XZ = imresize(flipud(permute(SI(:,round(size(SI,2)/2),:,:),[3 1 4 2])),scale); % Front
    XY = imresize(flipud(permute(SI(:,:,round(size(SI,3)/2),:),[2 1 4 3])),scale); % Top
end

function [XZ, YZ, XY]=magnitude_map(I)
    C=70/1000;
   % I is a 4D image {Z, X, Y, RS}
    
    SI = (sum(I(:,:,:,1:3).^2,4).^(1/2)).*C; % Magnitude matrix
    SI(I(:,:,:,4)<0)=Inf; % Non-convergent
    %SI(I(:,:,:,4)<-3)=Inf; % Masked
    
    scale=1;
    
    % extract slices
    YZ = imresize(flipud(permute(SI(round(size(SI,1)/2),:,:),[3 2 1])),scale); % Side
    XZ = imresize(flipud(permute(SI(:,round(size(SI,2)/2),:),[3 1 2])),scale); % Front
    XY = imresize(flipud(permute(SI(:,:,round(size(SI,3)/2)),[2 1 3])),scale); % Top
    
end

function [SI] = histo_Axis(I)

    C=70/1000;
   % I is a 4D image {Z, X, Y, RS}
    SI = I(:,:,:,1:3).*C; % Magnitude matrix    
    nc=(I(:,:,:,4)<0); % Non-convergent

    for c=1:3 % Loop over Components  
        sc=I(:,:,:,c);
        sc(nc)=nan;
        SI(:,:,:,c)=sc; % Non-convergent    
    end
end

function [] = plot_group(S,I,tit)

fig=figure();

%ColorMapped
for sp=1:3
    subplot(3,3,sp)
    imshow(S{sp})
end

%Magnitude
%Cmap= [ 0 1 0; flipud(hot(100)) ; 0 0 1];
Cmap= [ flipud(hot(20)) ; 0 0 1];

for sp=4:6
    subplot(3,3,sp)
    imagesc(S{sp})
    colormap(Cmap)
end
    colorbar

            
%Histogram
order=[1 3 2];
for sp=7:9
    m=I(:,:,:,order(sp-6));
    subplot(3,3,sp)
    histogram(m(~isnan(m)), 50, 'Normalization','PDF','BinMethod','scott');
    grid on;
    xlim(quantile(m(:),[0.05 0.95]));
end

sgtitle(tit)
%saveas(gcf,[tit '.png'])
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 14 12];
print(fig,tit,'-dpng')
print(fig,tit,'-dtiffn')

end