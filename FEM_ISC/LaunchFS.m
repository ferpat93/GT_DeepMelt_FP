% Launch Feature selection algorithm
%Load Data
warning('off','all')

%[Cases_filename,Cases_filepath] = uigetfile('*.csv','Select Cases table','C:\Users\lfp3\Dropbox\GT\Spring-18\IS_Paper\Matlab\CombinatorialCases.csv'); %PC
[Cases_filename,Cases_filepath] = uigetfile('*.csv','Select Cases table','/Users/lfp3/Dropbox/GT/Spring-18/IS_Paper/Matlab/CombinatorialCases.csv'); %MAC 


Parameters_Array = csvread(fullfile(Cases_filepath,Cases_filename),1,1);
nCases=length(Parameters_Array(:,1));
Parameters=[Parameters_Array(:,1)./1e7 Parameters_Array(:,2) tan(Parameters_Array(:,3)) tan(Parameters_Array(:,5))];

%Resp = load(fullfile('C:\Users\lfp3\Dropbox\GT\Spring-18\IS_Paper\Matlab',filesep,'Y_Responses.mat')); % PC
Resp = load(fullfile('/Users/lfp3/Dropbox/GT/Spring-18/IS_Paper/Matlab',filesep,'Y_Responses.mat')); % MAC
Cav_Responses=Resp.YCav;
EP_Responses=Resp.YEP;

Cav_Out=cell(4,2);
EP_Out=cell(4,2);

% Run for Cavity
disp('Cavity')
for i=1:length(Cav_Responses(1,:)) 
    disp(strcat('Column: ',num2str(i)));
    [Cav_Out{i,1}, Cav_Out{i,2}]=FeatSel(Parameters,Cav_Responses(:,i)); 
end

 %Run for EP
 disp('EP')
for i=1:length(EP_Responses(1,:))
    
    disp(strcat('Column: ',num2str(i)));
    [EP_Out{i,1}, EP_Out{i,2}]=FeatSel(Parameters,EP_Responses(:,i));   
end



zx