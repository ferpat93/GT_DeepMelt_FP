figure(10)
set(gcf,'units','normalized','outerposition',[0 0 1 1])

for t = 1:length(Treatments)
    hold on
    y=100-tot{1,1,t}(x)-tot{1,1,t}(x);
    plot(x,y);
end

legend(Treatments,'Position',[0.8 0.5 0.1 0.05])

%Save figures to desktop
saveas(gcf,fullfile(ResultsFolder,'StatisticAnalysis-Agar.tif'),'tiffn')