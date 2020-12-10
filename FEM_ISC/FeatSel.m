% Sequential feature selection

function [inmodel, history]=FeatSel(X,Y)
%% Define options for algorithm
opts = statset('display','iter','UseParallel',true);  

%% Call feature selection algorithm

[inmodel, history] = sequentialfs(@critfun,X,Y,...
                       'cv','none',...
                       'options',opts,...
                       'direction','forward',...
                       'nfeatures',4);
 
end

%% Function to evaluate accuracy

function dev = critfun(X,Y)
    %model=MultiPolyRegress(X,Y,2);
    %tbl = table(X(:,1),X(:,2),X(:,3),X(:,4),Y,'VariableNames',{'E','v','phi','psi','Y'});
    %model = fitlm(tbl,'Y ~ E*v*phi*psi - 1');
    model = fitlm(X,Y,'quadratic');
    assignin('base','model',model)
    dev = (model.MSE);
    %dev = 1-(model.Rsquared.Ordinary);
    
end
                   