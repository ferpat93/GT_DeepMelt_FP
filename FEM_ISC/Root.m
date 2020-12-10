 % Root Code for IS Paper 2020

%% 1) Set Cases and preprocessing

%%% INP File

RootFolder= pwd;
% [INP_name,INP_Path] = uigetfile('*.inp','Select the Base Input file','C:\Users\lfp3\Dropbox\GT\Spring-18\IS_Paper\Matlab\PS-Quadratic.inp');
% RootINP_fullpath=[INP_Path INP_name]; %fullpath of Initial INP file
RootINP_fullpath = '/Users/lfp3/Dropbox (GaTech)/GT/Fall_19/IS_2020/Matlab/Matlab\PS-Quadratic.inp';

%%% Parameter Bounds

%[Cases_filename,Cases_filepath] = uigetfile('*.csv','Select Cases table','C:\Users\lfp3\Dropbox\GT\Fall_19\IS_2020\Matlab\ParameterRangesAbaqus.csv');
Parameters_path = '/Users/lfp3/Dropbox (GaTech)/GT/Fall_19/IS_2020/Matlab/ParameterRangesAbaqus.csv';
Parameters_Array = csvread(Parameters_path,1,1);
Parameters_Table = readtable(Parameters_path);

ComputationPath='/Users/lfp3/Dropbox (GaTech)/GT/Fall_19/IS_2020/Computations';
%LinestoChangeINP=[7598,7594,7596,7626,7666];
LinestoChangeINP=[14130,14126,14128,14158,14204];

%% 2) Parameters

%Fixed Parameters
ratio=1; %??
yield=30000; %??
sp=500e3;

%% LOOP 
nSamples=1000; % per vertical stress
varNames = {'E','v','phi','SV','E1','E2','S1','S2'};
varTypes = {'string','double','double','double','double','double','double','double','double'};
Data = zeros(nSamples*Parameters_Array(4,3), numel(varNames));

VertStresses = sp.*linspace(Parameters_Array(4,1),Parameters_Array(4,2),Parameters_Array(4,3));

for svi =1:numel(VertStresses)
    SV = VertStresses(svi);
    parfor s=1:nSamples
        i = (svi-1)*nSamples + s;
        
        NameCase= strcat('C',num2str(i)); %Gets the name of the Case considered
        folderCase=fullfile(ComputationPath, NameCase); %Full path of the folder to be created

        % Run abaqus and extract the EP boundary and Cavity Shape
        Parameters_i = [(unifrnd(Parameters_Array(1:end-1,1),Parameters_Array(1:end-1,2)))' SV];
        Indexes_i=AbaqusRunExtract(RootINP_fullpath,NameCase,folderCase,Parameters_i,LinestoChangeINP); %Creates a modified input file for the case and returns its fullpath 
    end
end
