% Confusion

load('SVM')

X=table2array(SVM.ClassificationSVM.X);
Yt = SVM.ClassificationSVM.Y;
Yp = SVM.predictFcn(X);

t=[1 0 0; 0 1 0; 0 0 1];
YTT=t(Yt,:);
YPP=t(Yp,:);

plotroc(YTT',YPP')

%plotconfusion(YTT',YPP')

%set(gca,'xticklabel',{'Low' 'Mid' 'High' ''})
%set(gca,'yticklabel',{'Low' 'Mid' 'High' ''})