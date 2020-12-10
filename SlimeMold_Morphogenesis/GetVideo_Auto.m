%Create video
function []=GetVideo_Auto(Images,GR,rows,columns,Folder,Name)
warning off all

%Color Code
Colors=[206 185 28; 255 255 255 ; 183 219 127; 0 0 0 ].*(1/255);

%Video Parameters

filepathE = fullfile(Folder,strcat(Name,'-Video'));

FPS=4; % Frames per second of video
lengthV=90; %Total video length (seconds)
nImages=length(Images(1,:));
nF=min(nImages,round(FPS*lengthV)); %Total number of frames to display
sp=floor((nImages-1)/(nF-1)); %Spacing between displayed frames
%nFram=1+floor((nF-1)/sp);
F(1) = struct('cdata',[],'colormap',[]); % Creates struct variable to store figure frames

%% First Video - Entities

indexFrame=1;
for i=1:sp:nImages    
    I_list=zeros(rows*columns,3);    
    I_list(:,1:3)=Colors(Images(:,i),:);
    I_matrix=reshape(I_list,rows,columns,3);
    
    figure(17) 
    
    axis equal % Sets equal scale for axis

    image(I_matrix)
    drawnow %Update plot
    
    F(indexFrame) = getframe(gcf); %Get frame and store it    
    indexFrame=indexFrame+1;
end

assignin('base', 'I_matrix', I_matrix)

v = VideoWriter(filepathE); % Write video in the desired location
v.Quality=100;
v.FrameRate = round(FPS); %Sets frames per second
open(v) %Activate video
writeVideo(v,F) %Save it
close(v) %Deactivate video


%% Second Video - Entities

