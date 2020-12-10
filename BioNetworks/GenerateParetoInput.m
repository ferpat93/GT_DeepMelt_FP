function GenerateParetoInput(points,name)

path_terminal='/Users/lfp3/Dropbox\ \(GaTech\)/GT/Previous\ Semesters/Summer-17/BioNetworksPaper/src-website';
path_matlab='/Users/lfp3/Dropbox (GaTech)/GT/Previous Semesters/Summer-17/BioNetworksPaper/src-website';

DataFolder=fullfile(path_matlab,'data',name);
mkdir(DataFolder);

% Points [nX2] = list of coordinates X Y for points.
% First point corresponds to the root (initial point)
% All subsequent points are leaves.

%Stores text file to the location given by path


nPoints=size(points,1);

%% Create .cp.txt File

comma=',';
endline=',0.0';

cpPath=fullfile(DataFolder,[name '.cp.txt']);
fout = fopen(cpPath,'w');

line1=['"r -nom-",' num2str(points(1,1)) comma num2str(points(1,2)) endline]; % First Line, root
fprintf(fout,'%s\n',line1);

% Keep going with leaves 

for i=2:nPoints
    line=['"c' num2str(i-1) ' -nom-",' num2str(points(i,1)) comma num2str(points(i,2)) endline]; % First Line, root
    fprintf(fout,'%s\n',line);
end

% Create false connections

fprintf(fout,'%s\n',' ');
fprintf(fout,'%s\n','#edges');

s='r: ';

for i=2:nPoints
    s=[s 'c' num2str(i-1) ',']; % append line
end

fprintf(fout,'%s\n',s(1:end-1));

fclose(fout);

%% Create .sh file

OriginalSH=fullfile(path_matlab,'run_single.sh');
ModifiedSH=fullfile(path_matlab,'data',name,[name '.sh']);

% Start modifying text
fin = fopen(OriginalSH,'r');
fout = fopen(ModifiedSH,'w');

% Path for terminal
cpTerminal=fullfile('data',name,[name '.cp.txt']);

idk=0;
while ~feof(fin)
  
    idk=idk+1;
    s = fgetl(fin);
    
    if contains(s,'FILENAME="')
        s=['FILENAME="' cpTerminal '"'];
    end
    
    fprintf(fout,'%s\n',s);
end

fclose(fin);
fclose(fout);


%% Run Navlakha Script
%cd(DataFolder)

TerminalSH=['. ' path_terminal,'/data/',name,'/',[name '.sh']];
disp(TerminalSH);

clipboard('copy',TerminalSH)



[status,console] = system(TerminalSH)

if status~=0
    disp('There was an error')
end


end
