% Format: 
% [ Treatment Replicate Time GRE SM Ag R]    
load('E:\Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS\RESULTS\DecisionData.mat')

Treatments={'Glucose_100mM','Glucose_200mM','Control','NaCl_100mM'};
% GR Entities 
        % -1 -> Refinement
        % 0 -> No change
        % 1 -> Primary Growth
        % 2 -> Secondary Growth
GRES=1:3; % positives due to uint8 type
nS=9; % neighbourhood size

for treat=1:numel(Treatments)
    % [GRE SM Ag R] 
    Mt=masterArray((masterArray(:,1)==treat),4:7);    
    for GRE=1:numel(GRES)
        M=Mt((Mt(:,1)==GRES(GRE)),2:4);
        figure(1)
        hist=histogram(double(M(:,2))./sum(M(:,2:3),2),8,'Normalization','probability');
        y=hist.Values;
        x=hist.BinEdges(1:end-1)+(0.5*hist.BinWidth);
        figure(treat+100)
        hold on
        plot(x,y)
        title(Treatments(treat))
        legend('Refinement','Primary Growth','Secondary Growth')
    end       
end