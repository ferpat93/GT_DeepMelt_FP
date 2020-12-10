% trial
t=(2:420)*5/60;
Treatments={'Glucose_100mM','Glucose_200mM','Control','NaCl_100mM'};

for treat=1:4
    
    figure(45)
    subplot(2,2,treat)
    hold on
    plot(t,smooth(Data(2:end,1,treat)))
    plot(t,smooth(Data(2:end,4,treat)))
    legend('Choice','Random')
    title(Treatments{treat})
    
    figure(96+treat)
    hold on
    lt={'-',':',':'};
    c={'b','r'};
    ind=[1 4];
    for e=1:2
        for L=1:3
            plot(t,smooth(Data(2:end,ind(e)+L-1,treat)),strcat(lt{L},c{e}))
        end
    end
    
    legend('Choice','1st Quart Choice','3rd Quart Choice','Random','3rd Quart Random','3rd Quart Random')
    title(Treatments{treat})
    
    
    figure (10)
    hold on
    plot(t,Data(2:end,7,treat))
    
    figure(20)
    hold on
    plot(t,Data(2:end,1,treat)./Data(2:end,4,treat))

end

figure (10)
legend(Treatments)
title('Extent of growth region')

figure(20)
legend(Treatments)
title('Ratio of probabilities (>1 means unexplored is preferred')
