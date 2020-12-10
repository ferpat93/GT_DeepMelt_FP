close all;

DataFolder = '/Users/lfp3/Dropbox (GaTech)/GT/Spring-20/CE_paper/2ND-WEEK-PVM/';
%ResultsFolder = '/Users/lfp3/Dropbox (GaTech)/GT/Spring-20/SMN/Data';
k=-2; % Constant of pressure sensor.

Tests=GetSubfolders(DataFolder);
spl=cellfun(@(x) strsplit(x,'-'),Tests,'UniformOutput',false);
T = cell2table(vertcat(spl{:}),'VariableNames',{'Number','Length','Load','Density'});
T.Number = cell2mat(cellfun(@(x) str2double(x(2:end)),T.Number,'UniformOutput',0));
T.Length = cell2mat(cellfun(@(x) str2double(x(2:end)),T.Length,'UniformOutput',0));
T.Name = Tests;
T = sortrows(T,'Number');

CCfile = 'CC.txt';
CFfile = 'CF.txt';
SVfile = 'SV100.txt';

peaks = zeros(numel(Tests),2); % max load and vol at max
PV = cell(numel(Tests),1);
CC = cell(numel(Tests),1);
CF = cell(numel(Tests),1);

for f=1:size(T,1)
    try
        CC{f} = dlmread(fullfile(DataFolder,T.Name{f},CCfile));
        CF{f} = dlmread(fullfile(DataFolder,T.Name{f},CFfile));
    end
    
    SV = dlmread(fullfile(DataFolder,T.Name{f},SVfile));
    PV{f} = interpolatedPV(SV,CC{f},CF{f},k);
    
    [m, I] = max(PV{f}(:,2));
    peaks(f,:) = [m PV{f}(I,1)];

end

%% PV PLOTS
 
TestVariables = [ [2 4 6 2 4 6 2 4 2 6 6 4 6]' [1 1 1 1 1 1 0 0 0 0 0 0 0]' [0 0 0 1 1 1 0 0 1 0 1 1 2]']; % Length (2,4,6), density (0,1), load(0,1,2)

v=0:(100/58):100;
pPS = zeros(13,2);
testsToPlot = [9 12 11];
ImagesToPlot = cell(2,numel(testsToPlot),2);
maxStrain = [0 0]; steptoplot = 1;

colorlines = {'r','k','b'};
linetyp = {'-.','-'};
lims = [170 100 100];

f=figure;
for l=2:2:6 %loop lengths
    
    idx = find(TestVariables(:,1)==l);

    subplot(1,3,l/2); hold on;
    
    for i=1:numel(idx) % all tests of such length
        
        PVi = PV{idx(i)};
        plot(PVi(:,1),PVi(:,2),[colorlines{TestVariables(idx(i),3)+1} linetyp{TestVariables(idx(i),2)+1}])

    end
    ylim([0 lims(l/2)])
    grid on
    xlabel('Volume [cm^3]'); 
    ylabel('Pressure [kPa]'); 
    title(strcat('L =',{' '},num2str(l),'D')); grid on

end

figure(f); subplot(1,3,3); legend({'Dense - \sigma1','Dense - \sigma2','Loose - \sigma1','Loose - \sigma2','Loose - \sigma3'},'Location','southoutside');


%% Relationships Plot

stresses = [3.5 6 8]; densities = [0.4 0.75];
TV = [TestVariables peaks(:,1)];
TV(:,2) = densities(TestVariables(:,2)+1)';
TV(:,3) = stresses(TestVariables(:,3)+1)';


figure;
scatter(sv,peaks(:,1))


aaa = peaks(:,1)./sv;

%% CALIBRATION Example

test_to_plot = 1;
f=figure;
% CC - Constrained Calibration
subplot(2,1,1)
plot(CC{test_to_plot}(:,2),k.*CC{test_to_plot}(:,3),'k')
xlabel('Volume [cm^3]'); ylabel('Pressure [kPa]')
title('Constrained Expansion')
xlim([0 1.25]); ylim([0 210]); grid on

% CF - Free Calibration
subplot(2,1,2)
plot(CF{test_to_plot}(:,2),k.*CF{test_to_plot}(:,3),'k')          
xlabel('Volume [cm^3]'); ylabel('Pressure [kPa]')
title('Unconstrained Expansion')
xlim([0 2]); ylim([0 37.5]); grid on



%% Auxiliar Functions

function [PV] = interpolatedPV(SV,CCi,CFi,k)

    % Volume (1) - Pressure (2)
    [PV(:,1),i1] = unique(SV(:,2),'rows','stable');  PV(:,2) = k.*(SV(i1,3));   
    
    maxvol = min(max(PV(:,1)),max(CFi(:,2)));
    v = (0:0.005:maxvol)'; [~,im]=max(PV(:,1));
    p = interp1(PV(1:im,1),PV(1:im,2),v);
    
    if isempty(CFi)
        PV = [v p];
        return
    end
    val=find(CFi(:,2)>=max(v),1);
    [CF(:,1),i1] = unique(CFi(1:(val-1),2),'rows','stable');  CF(:,2) = k.*(CFi(i1,3));   
    pc = p - interp1(CF(:,1),CF(:,2),v); % Corrected pressure
    
%     val=find((CC(:,3).*k)>max(p),1);
%     CC = unique([CC(1:(val-1),2) k.*(CC(1:(val-1),3))],'rows','stable');

    PV = [v pc];
end





%% AUXILIAR FUNCTIONS

function [names]=GetSubfolders(path)
        d = dir(path);
        isub = [d(:).isdir]; %# returns logical vector
        names = {d(isub).name}';
        names(ismember(names,{'.','..'})) = [];
end 