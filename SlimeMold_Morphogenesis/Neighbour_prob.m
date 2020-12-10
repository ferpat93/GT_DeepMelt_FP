% Neighborhood probability

% Get a list of all files and folders in this folder.
Treatments={'Glucose_100mM','Glucose_200mM','Control','NaCl_100mM'};
root='E:\';
ResultsFolder=fullfile(root,'\Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS\RESULTS\');
files = dir(ResultsFolder);
% Extract only those that are directories.
subFolders = files([files.isdir]);
% Print folder names to command window.

% Loop over folders
Counter=zeros(4,1); % Treatment counter
time_range=[250 350];

masterArray = int8.empty(0,7);

for k = 3 : length(subFolders)
    %disp('Folder: ')
    k
    
    Name=subFolders(k).name(9:end);
    Specific_ResultsPath=fullfile(ResultsFolder,subFolders(k).name);
    
    for treat=1:4 % For loop to find to which treatment the folder belongs
        if contains(Name,Treatments{treat})
            Counter(treat)=Counter(treat)+1;
            rep=Counter(treat);    
            break
        end
    end
    
    [Entities,GR]=GetEntitiesGR(Specific_ResultsPath,Name,time_range); % Entities are time shifted already

    % GR Entities 
        % -1 -> Refinement
        % 0 -> No change
        % 1 -> Primary Growth
        % 2 -> Secondary Growth
  
     nT=size(GR,3);
     
     % Format: 
     % [ Treatment Replicate Time GRE SM Ag R]
     nEntries=sum(abs(GR)>0,[1 2]); % Points of interest per time step
     masterArray_f=ones(sum(nEntries),7,'int8');
     masterArray_f(:,1)=treat; % Asign Treatment
     masterArray_f(:,2)=rep; % assing Replicate

    initial_row=1;
    for t=1:nT
        end_row=initial_row+nEntries(t)-1;
        masterArray_f(initial_row:end_row,3)=t; % Asign Treatment
        masterArray_f(initial_row:end_row,4:7)=getNeigh(GR(:,:,t),Entities(:,:,t));   
        initial_row=end_row+1;
    end
     
    masterArray=[masterArray;masterArray_f];
   
end



 
%% Auxiliary Functions

function [ImagesArray, GR]=GetEntitiesGR(path, name,time_range)
    
    images_path =fullfile(path,strcat(name,'-ClassifiedImages'));
    struct = load(images_path);   
    ImagesArray=struct.ImagesArray;
    rows=struct.rows;
    columns=struct.columns;
       
    if numel(unique(ImagesArray(:,200)))<4
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
    
    GR_path =fullfile(path,strcat(name,'-GR_Images'));
    
    try 
        GRstruct = load(GR_path);
        GR=GRstruct.GR;
    catch 
        %Growth-Refining-Superpose Method
        GR=zeros(size(ImagesArray));
        aa=ImagesArray.^2;
        delta=aa-[aa(:,1) aa(:,1:end-1)];
        GR(delta==8)=-1; % Refining
        GR(delta==-8)=2; % Secondary growth (superpose)
        GR(delta==-3)=1; % Primary Growth
        
        save(GR_path,'GR','-v7.3');
    end
   
    % Resize Image and GR arrays
    ImagesArray=reshape(ImagesArray(:,time_range(1)-1:time_range(2)-1),[rows columns (time_range(2)-time_range(1)+1)]);
    GR=reshape(GR(:,time_range(1):time_range(2)),[rows columns (time_range(2)-time_range(1)+1)]);
 
end

function [M]=getNeigh(GR,Entities)
    M=zeros(sum(abs(GR)>0,'all'),4,'int8');
    ln=2;%padding around point
    
    % GR Entities 
        % -1 -> Refinement
        % 0 -> No change
        % 1 -> Primary Growth
        % 2 -> Secondary Growth
    GRES=[-1,1,2]; 
    
    % Trinary Entities 
        % 1 -> Slime Mold
        % 2 -> Agar
        % 3 -> Residumm
        % 4 -> Mask (outside dish)
        
    ENTS=1:3;  
    in_row=1; 
    for GR_i=1:numel(GRES) % For Every change entity
        [rows, columns]=find(GR==GRES(GR_i)); % Find the points of occurence
        M_i=zeros(numel(rows),4,'int8');
        M_i(:,1)=GR_i;
        for p=1:numel(rows) % for every point find the neighbourhood composition
            Neigh=Entities(rows(p)-ln:rows(p)+ln,columns(p)-ln:columns(p)+ln);
            M_i(p,2:end)=[sum(Neigh==ENTS(1),'all') sum(Neigh==ENTS(2),'all') sum(Neigh==ENTS(3),'all')];  
        end
        end_row=in_row+numel(rows)-1;
        M(in_row:end_row,:)=M_i;
        in_row=end_row+1;
    end
end















