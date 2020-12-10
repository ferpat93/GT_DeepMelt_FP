function []=AbaqusRunExtract(RootINP_fullpath,NameCase,FolderName,parameters,index)

% index - Named LinestoChangeINP, index of lines to change in the text
% (INP) file

mkdir(FolderName) %creates folder for the file
cd(FolderName)

%% MODIFY INP FILE
Case_INP=strcat(NameCase,'.inp');

% Start modifying text
fin = fopen(RootINP_fullpath,'r');
fout = fopen(Case_INP,'w');

idk=0;
while ~feof(fin)
    
    idk=idk+1;
    s = fgetl(fin);
    
    if and(idk>=min(index),idk<=max(index))
        switch idk
            case index(1) % Elastic Parameters
                s=strcat(num2str(1000.*parameters(1)),'., ',num2str(parameters(2)));
            case index(2) % Friction/dilatancy angles - DP
                s=strcat(num2str(parameters(3)),'., ',num2str(parameters(4)),'., ',num2str(parameters(5)));
            case index(3) % Yield Stress
                s=strcat(num2str(parameters(6)),',0.');
            case index(4) % Vertical Stress
                s=strcat('Top_Surface, P, ',num2str(parameters(7)));
            case index(5) % Cavity Stress
                s=strcat('CavityWall_Surface, P, ',num2str(parameters(8)));         
        end
    end
fprintf(fout,'%s\n',s); 
end

fclose(fin);
fclose(fout);

%copyfile([Inp_file '.inp'],[pwd '\' Inp_file '\'])

         
%% RUN ABAQUS - FROM INP 
      
disp('Simulation Started') 
tic
% Run the input file with Abaqus
system(['abaqus job=' NameCase]);

%%%%abaqus viewer database=cantilever

sw=true; 
tic; 
while sw 
    % Pause Matlab execution in order for the lck file to be created 
    pause(0.5); 
    % While the lck file exists, pause Matlab execution. If it is 
    % deleted, exit the while loop and proceed. 
    while exist([NameCase '.lck'],'file')==2 
        pause(0.1) 
        % the lck file has been created and Matlab halts in this loop. 
        % Set sw to false to break the outer while loop and continue 
        % the code execution. 
        sw=false; 
    end 
    % In case that the lck file cannot be detected, then terminate 
    % infinite execution of the outer while loop after a certain 
    % execution time limit (5 sec) 
    if sw && (toc>5) 
        sw=false; 
    end 
end 
% NOTE: Alternatively, you can replace lines 27 to 49 by system(['abaqus job=' Inp_file ' interactive']) 
 
disp('Simulation Finished') 
toc

%% RETRIEVE CONTOUR MAP

% Run the Python file from the shell

system(['abaqus cae noGUI==' Export.py]);
