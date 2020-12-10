%% Generate INP Files

RootFolder= pwd;
[INP_name,INP_Path] = uigetfile('*.inp','Select the Base Input file','C:\Users\lfp3\Dropbox\GT\Spring-18\IS_Paper\Matlab\PS-Quadratic.inp');
RootINP_fullpath=[INP_Path INP_name]; %fullpath of Initial INP file
[Cases_filename,Cases_filepath] = uigetfile('*.csv','Select Cases table','C:\Users\lfp3\Dropbox\GT\Spring-18\IS_Paper\Matlab\CombinatorialCases.csv');
Parameters_Array = csvread(fullfile(Cases_filepath,Cases_filename),1,1);
Parameters_Table = readtable(fullfile(Cases_filepath,Cases_filename));
nCases=length(Parameters_Array(:,1));

ComputationPath='C:\Users\lfp3\Documents\IS_SensitivityAnalysis\Cases\INP\';
LinestoChangeINP=[14130,14126,14128,14158,14204];

cd(ComputationPath)
%% LOOP 
for Case=1:nCases
    
    % Import Parameters
    NameCase=Parameters_Table.Parameter{Case}; %Gets the name of the Case considered
    Parameters_i=Parameters_Array(Case,:);

    %% MODIFY INP FILE
    Case_INP=strcat(NameCase,'.inp');

    % Start modifying text
    fin = fopen(RootINP_fullpath,'r');
    fout = fopen(Case_INP,'w');

    idk=0;
    while ~feof(fin)
        idk=idk+1;
        s = fgetl(fin);
        if and(idk>=min(LinestoChangeINP),idk<=max(LinestoChangeINP))
            switch idk
                case LinestoChangeINP(1) % Elastic Parameters
                    s=strcat(num2str(Parameters_i(1)),' , ',num2str(Parameters_i(2)));
                case LinestoChangeINP(2) % Friction/dilatancy angles - DP
                    s=strcat(num2str(Parameters_i(3)),' , ',num2str(Parameters_i(4)),' , ',num2str(Parameters_i(5)));
                case LinestoChangeINP(3) % Yield Stress
                    s=strcat(num2str(Parameters_i(6)),',0.');
                case LinestoChangeINP(4) % Vertical Stress
                    s=strcat('Top_Surface, P, ',num2str(Parameters_i(7)));
                case LinestoChangeINP(5) % Cavity Stress
                    s=strcat('CavityWall_Surface, P, ',num2str(Parameters_i(8)));         
            end
        end
    fprintf(fout,'%s\n',s); 
    end

    fclose(fin);
    fclose(fout); 
    
end
