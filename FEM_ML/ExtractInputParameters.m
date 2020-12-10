% Extract parameters from inp file %
%{
*Material, name=Sand
*Density
 2276.53,
*Drucker Prager
 38.7997,      1., 15.5145
*Drucker Prager Hardening, type=SHEAR
 7468.7,0.
*Elastic
 5.08502e+07, 0.388685
 
*Initial Conditions, type=STRESS, GEOSTATIC
%}

Lines_to_find = {'*Density','*Drucker Prager','*Drucker Prager Hardening',...
                '*Elastic','*Initial Conditions,'};
            
path = pwd;
path = 'C:\Users\lfp3\Documents\Local_AutoAbaqus';
Type = 'Disp_'; %'Pres_';

Folders = GetSubfolders(path,Type);

ParArray = zeros(numel(Folders),7);
prefix_inp = '_Standard.inp';

for f=1:numel(Folders)
    
    %% MODIFY INP FILE
    INP_file=fullfile(path,Folders{f},strcat(Folders{f},prefix_inp)); % INP path
    fin = fopen(INP_file,'r'); % Open INP

    idk=0; % current line being read
    search_line = 1; % Current line being searched for
    
    GetLine = false;
    while and(~feof(fin),search_line<=numel(Lines_to_find))

        idk=idk+1; % Next Line
        s = fgetl(fin); % Read line
        
        
        if GetLine
            l=str2double(strsplit(s,',')); % Parse line
            switch search_line
                case 1 % Density
                    d = l(1);
                case 2 % Friction/dilatancy angles - DP 
                    psi = l(3);
                    beta = l(1);
                case 3 % Yield Stress - DP
                    y = l(1);
                case 4 % Elastic (E,v)
                    E = l(1); v = l(2);
                case 5 % Depth (H)
                    H=abs(l(5))-15;                
            end
            search_line = search_line + 1; % Next line
        end

        GetLine = false;
        if search_line<=numel(Lines_to_find)
            if contains(s,Lines_to_find{search_line})
                GetLine = true;
            end
        end
    end

    fclose(fin); % Close text file
    
    % Translate parameters
    tanpsi = tand(psi);
    M = (3*(9-tanpsi^2))^(1/2)/(9-tanpsi*tand(beta));
    phi = asind(M*tand(beta));
    c = M*y/cosd(phi) ;

    ParArray(f,:) = [H d E v phi psi c]; % Store parameters
    
end

TableHeaders = {'Depth','d','E','v','phi','psi','c'};
Table = array2table(ParArray);
Table.Properties.VariableNames = TableHeaders;
Table.Folder = Folders;
Table = [Table(:,end) Table(:,1:end-1)];

writetable(Table,[Type 'InputParameters.txt']);


%% Auxiliar Functions

function [names]=GetSubfolders(path,string)
        d = dir(path);
        isub = [d(:).isdir]; %# returns logical vector
        names = {d(isub).name}';
        %names(or(ismember(names,{'.','..'}),~contains(names,string))) = [];
        names(~contains(names,string)) = [];
end 