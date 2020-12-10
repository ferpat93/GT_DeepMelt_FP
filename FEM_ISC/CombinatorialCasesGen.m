% Generate Combinatorial Cases

%Import the file with the ranges

[Cases_filename,Cases_filepath] = uigetfile('*.csv','Select Cases table','C:\Users\lfp3\Dropbox\GT\Fall_19\IS_2020\Matlab\ParameterRangesAbaqus.csv');

Parameters_Array = csvread(fullfile(Cases_filepath,Cases_filename),1,1);
Parameters_Table = readtable(fullfile(Cases_filepath,Cases_filename));

ParLists=cell(length(Parameters_Array(:,1)),1);

%Fixed Parameters
ratio=1;
yield=30000;
sv=200e3;
sp=500e3;

nCases=1;
for i=1:length(ParLists)
     nCases=nCases*Parameters_Array(i,3);
end

E=linspace(Parameters_Array(1,1),Parameters_Array(1,2),Parameters_Array(1,3));
v=linspace(Parameters_Array(2,1),Parameters_Array(2,2),Parameters_Array(2,3));
teta=linspace(Parameters_Array(3,1),Parameters_Array(3,2),Parameters_Array(3,3));
Rpsi=linspace(Parameters_Array(4,1),Parameters_Array(4,2),Parameters_Array(4,3));

CombMatrix=cell(nCases,9);
Cases=zeros(nCases,4);
iC=0;
for iE=1:Parameters_Array(1,3)
    for iv=1:Parameters_Array(2,3)
        for it=1:Parameters_Array(3,3)
            for ip=1:Parameters_Array(4,3)
                iC=iC+1;
                %a=(Rpsi(ip))*(1+0.333*(sind(teta(it)))^2)^(-0.5);
                %b= atand(3^(0.5)*sind(teta(it)));
                b=teta(it);
                psi=b*Rpsi(ip);
                
                Cases(iC,1)=E(iE);  
                Cases(iC,2)= v(iv);   
                Cases(iC,3)=b;
                Cases(iC,4)=psi;   
                Cases(iC,5)=Rpsi(ip);   
                 
%                 CombMatrix{iC,1}=strcat('C',num2str(iC));   
%                 CombMatrix{iC,2}=E(iE);   
%                 CombMatrix{iC,3}=v(iv);   
%                 CombMatrix{iC,4}=b;   
%                 CombMatrix{iC,5}=ratio;  
%                 CombMatrix{iC,6}=psi;   
%                 CombMatrix{iC,7}=yield;   
%                 CombMatrix{iC,8}=sv;   
%                 CombMatrix{iC,9}=sp;   
                
            end        
        end        
    end    
end

gfds
%T = cell2table(CombMatrix,'VariableNames',{'Parameter','E','V','phi','ratio','psi','yield','sv','sp'});
T = cell2table(CombMatrix,'VariableNames',{'Parameter','E','V','phi','ratio','psi','yield','sv','sp'});
writetable(T,'CombinatorialCasesMC.csv');