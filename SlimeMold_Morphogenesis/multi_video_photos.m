%function []=multi_video_photos(f)

paths=cell(4,1);
srcFiles=cell(4,1);
imgType = '*.tif'; % change based on image type
maxRows=0;

for i=1:4
    B = strsplit(f{i},'\');
    paths{i}=fullfile(B{1},B{2},B{3},'PHOTOS',B{5}(9:end));
    srcFiles{i} =dir([paths{i} '/' imgType]); % Gets info of the files inside the given folder with the given extension
    filename1 = fullfile(paths{i},'/',srcFiles{i}(1).name); %get first image in the folder
    I = imread(filename1); % open the image
    maxRows=max(maxRows,size(I,1));
end


warning off all

% TEXT 
%text={'Glucose 100mM','Glucose 200mM','Control','NaCl 100mM'};
text={'Control + Glucose 100mM','Control + Glucose 200mM','NaCl 200mM + Glucose 100mM','NaCl 200mM + Glucose 200mM'};

x=205;
y=30;
position = [x+20 y; x+maxRows+20 y; x+10 y+hl+maxRows; x+maxRows-40 y+hl+maxRows]; 
box_color = {'white','white','white','white'};

%Video Parameters
filepathE = fullfile(pwd,'Spot-Photos-Video');

FPS=4; % Frames per second of video
lengthV=90; %Total video length (seconds)
nImages=420;
nF=min(nImages,round(FPS*lengthV)); %Total number of frames to display
sp=floor((nImages-1)/(nF-1)); %Spacing between displayed frames
%nFram=1+floor((nF-1)/sp);
F(1) = struct('cdata',[],'colormap',[]); % Creates struct variable to store figure frames

hl=150;
%TextRow=4.*ones(hl,2*maxRows,420);

%Images=uint8(zeros(2*maxRows+1+2*hl,2*maxRows+1,3,420));
row=[hl+1,     hl+1,      (2*hl)+maxRows+1, (2*hl)+maxRows+1];
col=[1,     maxRows+2,          1,              maxRows+2];

indexFrame=1;

for t=1:420
    
    Im=uint8(zeros(2*(maxRows+hl),2*maxRows+1,3));
    
    for i=1:4
        I = imread(fullfile(paths{i},'/',srcFiles{i}(t).name)); % open the image
        rows=size(I,1);
        I=uint8([I zeros(rows,maxRows-rows,3); zeros(maxRows-rows,maxRows,3)]);  
        Im(row(i):row(i)+maxRows-1,col(i):col(i)+maxRows-1,:)=I;
    end
    
    Im = insertText(Im,position,text,'FontSize',60,'BoxColor',box_color,'BoxOpacity',1);
    
    figure(17)   
    
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis equal % Sets equal scale for axis
    
    imshow(Im)
    drawnow %Update plot
    
    F(indexFrame) = getframe(gcf); %Get frame and store it    
    indexFrame=indexFrame+1;
    
end

v = VideoWriter(filepathE); % Write video in the desired location
v.Quality=100;
v.FrameRate = round(FPS); %Sets frames per second
open(v) %Activate video
writeVideo(v,F) %Save it
close(v) %Deactivate video


%% First Video - Entities
%rows=2*(maxRows+hl);
%columns=2*maxRows;
%{
indexFrame=1;
for i=1:sp:2%nImages    

    figure(17)   
    
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    axis equal % Sets equal scale for axis
    
    imshow(Images(:,:,:,i))
    drawnow %Update plot
    
    F(indexFrame) = getframe(gcf); %Get frame and store it    
    indexFrame=indexFrame+1;
end
%}
%{
%assignin('base', 'I_matrix', I_matrix)

v = VideoWriter(filepathE); % Write video in the desired location
v.Quality=100;
v.FrameRate = round(FPS); %Sets frames per second
open(v) %Activate video
writeVideo(v,F) %Save it
close(v) %Deactivate video
%}