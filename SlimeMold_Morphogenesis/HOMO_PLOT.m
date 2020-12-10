% FINAL Plot - SPOT ANALYSIS

% Check Progress of analyses
clear 
close all
clc

Treatments={'Glucose_100mM','Glucose_200mM','Control','NaCl_100mM'};
ResultsFolder='C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS\RESULTS\';

load(fullfile(ResultsFolder,'OutputData_HomoAnalysis')); % LOADS OUTPUT FILE
 
time=((0:5:449*5)/60)';
%% PLOT 1

%figure (1); % Proportion plots
 
set(gcf, 'Position', get(0, 'ScreenSize')); 
set(gcf,'name','SLIME MOLD - Image Analysis') 

colors={'r','b','k'};

for t=1:4
    %subplot(2,2,t)
    figure(10+t)
    title(Treatments{t})
     
    %set(gca, 'ColorOrder',lines(3),'NextPlot', 'replacechildren');
    
    hold on
    for ent=1:3       
        Data=reshape(Output{t}(ent,1:450,:),[450,length(Output{t}(ent,1,:))]);
        h=plot(time,nanmean(Data,2));
        set(h, {'color'}, colors(ent));
    end
    
    
    legend('Slime Mold','Agar','Residumm') 
    ylim([0 100])
    xlim([0 36]) 
    xlabel('Time of growth (Hours)')
    ylabel('Percentage of total Area')
    grid on
    grid minor
    
    
end

%saveas(gcf,fullfile(ResultsFolder,'ENTITIES_OVER_TIME.tif'),'tiffn')
%saveas(gcf,fullfile(ResultsFolder,'ENTITIES_OVER_TIME.png'),'png')


%% PLOT 2

%figure (2); % Growing-Refining fluctuations
 
set(gcf, 'Position', get(0, 'ScreenSize')); 
set(gcf,'name','Cumulative fluctuations') 


for t=1:4
    %subplot(2,2,t)
    figure(20+t)
    title(Treatments{t})
    hold on
    
    for l=1:3
        Data=reshape(Output{t}(3+l,1:450,:),[450,length(Output{t}(ent,1,:))]);
        D=cumsum(Data,1);
        h=plot(time,nanmean(D,2));
        set(h, {'color'}, colors(l));
    end
    
    legend({'Secondary Growth','Primary Growth','Refining'},'Location','northwest');
   
    xlim([0 40])
    ylim([0 600])
    xlabel('Time of growth (Hours)')
    ylabel('Cumulative percentage')
    grid on
    grid minor
    
end

saveas(gcf,fullfile(ResultsFolder,'GROWTH_FLUCTUATIONS.tif'),'tiffn')
saveas(gcf,fullfile(ResultsFolder,'GROWTH_FLUCTUATIONS.png'),'png')


%% PLOT 4

figure (4); % Shape plots
set(gcf, 'Position', get(0, 'ScreenSize')); % Maximize the figure. 
set(gcf,'name','Shape Indexes vs time') 

indexes={'Eccentricity','Solidity','Circularity','# Clusters'};

for index=1:4
     if index==3
        c=-1;        
    else
        c=1;       
    end
    for t=1:4
        subplot(4,4,(index-1)*4+t)
        title(Treatments{t});
        
        Data=reshape(Output{t}(6+index,1:450,:),[450,length(Output{t}(ent,1,:))]);

        hold on

        for s=1:length(Data(1,:))
            plot(time,Data(:,s).^c,':')
        end

        plot(time,nanmean(Data.^c,2))

        xlabel('Time of growth (Hours)')
        ylabel(indexes{index})
        xlim([0 40])
        %ylim([0 1])
    end
end

saveas(gcf,fullfile(ResultsFolder,'JOINT_SHAPE_INDEXES.tif'),'tiffn')
saveas(gcf,fullfile(ResultsFolder,'JOINT_SHAPE_INDEXES.png'),'png')
%% PLOT 5-8

for fig=5:8
    figure (fig); % Spatial plots
    % Maximize the figure. 
    set(gcf, 'Position', get(0, 'ScreenSize')); 
    set(gcf,'name',indexes{fig-4}) 
    
    if fig==7
        c=-1;        
    else
        c=1;       
    end
            
    for t=1:4
        subplot(2,2,t)
        title(Treatments{t});
        Data=reshape(Output{t}(fig+1,1:450,:),[450,length(Output{t}(ent,1,:))]);

        hold on

        for s=1:length(Data(1,:))
            
            plot(Data(:,s).^c,':')
        end

        plot(nanmean(Data.^c,2))

        xlabel('Frame')
        ylabel(indexes{fig-4})
        xlim([0 450])
        grid on
    end
    
    saveas(gcf,fullfile(ResultsFolder,strcat(indexes{fig-4},'.tif')),'tiffn')
end

%% Joined II

figure (10); % Shape plots
set(gcf, 'Position', get(0, 'ScreenSize')); % Maximize the figure. 
set(gcf,'name','Shape Indexes vs time') 

indexes={'Eccentricity','Solidity','Circularity','# Clusters'};

for index=1:4
        subplot(2,2,index)
       
        hold on
        for t=1:4
        Data=reshape(Output{t}(6+index,1:450,:),[450,length(Output{t}(ent,1,:))]);

        plot(time,nanmean(Data.^c,2))
        end
        xlabel('Time of growth (Hours)')
        ylabel(indexes{index})
        xlim([0 40])
        %ylim([0 1])

end

legend(Treatments,'location','best')

saveas(gcf,fullfile(ResultsFolder,'JOINED_SHAPE_INDEXES_II.tif'),'tiffn')
saveas(gcf,fullfile(ResultsFolder,'JOINED_SHAPE_INDEXES_II.png'),'png')