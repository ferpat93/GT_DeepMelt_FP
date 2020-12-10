
% Get Joint Video

f{2}='E:\Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS\RESULTS\Results-Glucose_200mM_R2_07_02_Et1\Glucose_200mM_R2_07_02_Et1-ClassifiedImages.mat';
f{3}='E:\Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS\RESULTS\Results-Control_R2_07_02_Et1\Control_R2_07_02_Et1-ClassifiedImages.mat';
f{1}='E:\Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS\RESULTS\Results-Glucose_100mM_R1_07_02_Et1\Glucose_100mM_R1_07_02_Et1-ClassifiedImages.mat';
f{4}='E:\Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS\RESULTS\Results-NaCl_100mM_R2_07_02_Et1\NaCl_100mM_R2_07_02_Et1-ClassifiedImages.mat';

I=cell(4,1);

maxRows=0;

for i=1:4
    struct = load(f{i},'rows');
    maxRows=max(maxRows,struct.rows);  
end

for i=1:4
    [A, rows, ~]=check_Trin(f{i});
    I{i}=[A 4.*ones(rows,maxRows-rows,420); 4.*ones(maxRows-rows,maxRows,420)]; 
end

hl=100;
TextRow=4.*ones(hl,2*maxRows,420);
BigImage=reshape([TextRow;I{1} I{2};TextRow;I{3} I{4}],[4*maxRows*(maxRows+hl) 420]);  

clear I A struct

warning off all
%Color Code
Colors=[206 185 28; 255 255 255 ; 183 219 127; 0 0 0 ].*(1/255);

% TEXT 

x=350;
y=15;
text={'Glucose 100mM','Glucose 200mM','Control','NaCl 100mM'};
position = [x y; x+maxRows y; x+150 y+50+maxRows; x+maxRows y+50+maxRows]; 
box_color = {'white','white','white','white'};

%Video Parameters
filepathE = fullfile(pwd,'Homo-Joint-Video');

FPS=4; % Frames per second of video
lengthV=90; %Total video length (seconds)
nImages=length(BigImage(1,:));
nF=min(nImages,round(FPS*lengthV)); %Total number of frames to display
sp=floor((nImages-1)/(nF-1)); %Spacing between displayed frames
%nFram=1+floor((nF-1)/sp);
F(1) = struct('cdata',[],'colormap',[]); % Creates struct variable to store figure frames

%% First Video - Entities
rows=2*(maxRows+hl);
columns=2*maxRows;

indexFrame=1;
for i=1:sp:nImages    
    I_list=zeros(rows*columns,3);    
    I_list(:,1:3)=Colors(BigImage(:,i),:);
    I_matrix=reshape(I_list,rows,columns,3);
    I_matrix = insertText(I_matrix,position,text,'FontSize',60,'BoxColor',box_color,'BoxOpacity',1);
    
    figure(17)   
    
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis equal % Sets equal scale for axis
    
    imshow(I_matrix)
    drawnow %Update plot
    
    F(indexFrame) = getframe(gcf); %Get frame and store it    
    indexFrame=indexFrame+1;
end

%assignin('base', 'I_matrix', I_matrix)

v = VideoWriter(filepathE); % Write video in the desired location
v.Quality=100;
v.FrameRate = round(FPS); %Sets frames per second
open(v) %Activate video
writeVideo(v,F) %Save it
close(v) %Deactivate video


%% Auxiliary Functions

function [ImagesArray, rows, columns]=check_Trin(path)
    
    struct = load(path);   
    ImagesArray=struct.ImagesArray;
    rows=struct.rows;
    columns=struct.columns;
    
    c=unique(ImagesArray(:,400));
    
    if numel(c)<4
        for p=1:rows*columns %Loop through every pixel
            firstC=find(ImagesArray(p,:)==1,1,'first');
            if or(firstC==columns,isempty(firstC)) % If there is never SM or only grows at the last picture
              continue
            else
              cs=find(ImagesArray(p,firstC+1:end)==2)+firstC; %find columns to change (equal to 2 after first SM)
              ImagesArray(p,cs)=3;
            end
        end
        save(path,'ImagesArray','rows','columns','-v7.3');
    end
    
    ImagesArray=ImagesArray(:,1:420);  
    ImagesArray=reshape(ImagesArray,[rows columns 420]);
  
end