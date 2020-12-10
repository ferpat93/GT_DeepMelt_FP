% SURVIVAL PLOT - SPOT ANALYSIS

Treatments={'Glucose_100mM_NaCl_200mM','Glucose_200mM_NaCl_200mM','Glucose_100mM','Glucose_200mM'};
ResultsFolder='C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS\';
load(fullfile(ResultsFolder,'OutputData_SpotAnalysis')); % LOADS OUTPUT FILE
 
%% Survival Plot 
filename=fullfile(ResultsFolder,'SurvivalData');

survivaldata=cell(4,1);
xlRange = 'A2';
tv=10:0.25:35;
PL=0;

for t=1:4
    Data=reshape(Output{t}(7,1:400,:),[400,20]);
    [fr, ~]= find(Data==0);   
    [p,ti,~,~] = ecdf(fr.*(5/60),'function','survivor');

%% WRITE EXCEL
    ti(2)=ti(2)-0.000000001;
    ff = interp1(ti,p,tv,'linear');%,'spline');
    survivaldata{t}=ff';
    xlswrite(filename,[tv' ff'],t,xlRange) 

    if PL==1
        figure(t)
        hold on
        plot(ti,p)    
        plot(tv,ff)
    end
end

sum(survivaldata{4}-survivaldata{3}>0)/sum(~isnan(survivaldata{4}-survivaldata{3}))
sum(survivaldata{2}-survivaldata{1}>0)/sum(~isnan(survivaldata{2}-survivaldata{1}))

sum(survivaldata{1}-survivaldata{3}>0)/sum(~isnan(survivaldata{1}-survivaldata{3}))
sum(survivaldata{2}-survivaldata{4}>0)/sum(~isnan(survivaldata{2}-survivaldata{4}))



