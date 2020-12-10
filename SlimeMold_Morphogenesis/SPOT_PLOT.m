% FINAL Plot - SPOT ANALYSIS

% Check Progress of analyses
clear 
close all
clc

Treatments={'Glucose_100mM_NaCl_200mM','Glucose_200mM_NaCl_200mM','Glucose_100mM','Glucose_200mM'};
ResultsFolder='C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS\';

load(fullfile(ResultsFolder,'OutputData_SpotAnalysis')); % LOADS OUTPUT FILE
 

%% PLOT 1

figure (1); % Proportion plots
 
set(gcf, 'Position', get(0, 'ScreenSize')); 
set(gcf,'name','SLIME MOLD - Image Analysis') 


for t=1:4
    subplot(2,2,t)
    title(Treatments{t})
      
    hold on
    for ent=1:3
        Data=reshape(Output{t}(ent,1:450,:),[450,20]);
        plot(nanmean(Data,2))
    end
    
    legend('Slime Mold','Agar','Residumm') 
    ylim([0 100])
    xlim([0 450]) 
    xlabel('Frame')
    ylabel('Percentage of total Area')
    grid on
    grid minor
end

saveas(gcf,fullfile(ResultsFolder,'ENTITIES_OVER_TIME.tif'),'tiffn')
%% PLOT 2

figure (2); % Growing-Refining fluctuations
 
set(gcf, 'Position', get(0, 'ScreenSize')); 
set(gcf,'name','Cumulative fluctuations') 


for t=1:4
    subplot(2,2,t)
    title(Treatments{t})
    hold on
    
    for l=1:3
        Data=reshape(Output{t}(3+l,1:450,:),[450,20]);
        D=cumsum(Data,1);
        plot(nanmean(D,2))
    end
    
    legend({'Primary Growth','Secondary Growth','Refining'},'Location','northwest');
    xlabel('Frame')
    ylabel('Cumulative percentage')
    xlim([0 450])
    ylim([0 425])
    grid on
    grid minor
end

saveas(gcf,fullfile(ResultsFolder,'GROWTH_FLUCTUATIONS.tif'),'tiffn')

%% PLOT 3

figure (3); % Spatial plots
% Maximize the figure. 
set(gcf, 'Position', get(0, 'ScreenSize')); 
set(gcf,'name','Distance to food vs time') 

for t=1:4
    Points=[];
    subplot(2,2,t)
    title(Treatments{t});
    Data=reshape(Output{t}(7,1:450,:),[450,20]);
    
    hold on
    
    for s=1:length(Data(1,:))
        plot(Data(:,s),':')
    end
    
    plot(nanmean(Data,2))
    
    xlabel('Frame')
    ylabel('Distance (pixels)')
    ylim([0 450])
    xlim([0 400])
    grid on
end

saveas(gcf,fullfile(ResultsFolder,'DISTANCE_TO_FOOD.tif'),'tiffn')


%% Survival Plot 
filename=fullfile(ResultsFolder,'SurvivalData');
colores=lines(4);
 
figure (10); % Spatial plots
% Maximize the figure. 
set(gcf, 'Position', get(0, 'ScreenSize')); 
set(gcf,'name','Time to reach Glucose') 

for t=1:4
    Points=[];
    subplot(2,2,t)
    
    Data=reshape(Output{t}(7,1:400,:),[400,20]);
    [fr, ~]= find(Data==0);
    
    [f,x,flow,fup] = ecdf(fr.*(5/60),'function','survivor');
    plot(x,f,'Color',colores(t,:));
    hold on
    plot(x,flow,':','Color',colores(t,:))
    plot(x,fup,':','Color',colores(t,:))

    xlabel('Growth time (Hours)')
    ylabel('P')
    legend(Treatments{t});
    ylim([0 1])
    grid on
end

saveas(gcf,fullfile(ResultsFolder,'SurvivalGrid.tif'),'tiffn')

figure (11); % Spatial plots
% Maximize the figure. 
set(gcf, 'Position', get(0, 'ScreenSize')); 
set(gcf,'name','Time to reach Glucose') 

%
survivaldata=cell(4,1);
xlRange = 'A2';
tv=10:0.25:35;
%

hold on

for t=1:4
    ct=colores(t,:);
    Points=[];
    Data=reshape(Output{t}(7,1:400,:),[400,20]);
    [fr, ~]= find(Data==0);
    
    [f,x,flow,fup] = ecdf(fr.*(5/60),'function','survivor');
    
    %% WRITE EXCEL
    survivaldata{t}=[x f];
    x(2)=x(2)-0.001;
    ff = interp1(x,f,tv,'spline');
    xlswrite(filename,[tv' ff'],t,xlRange) 
    %%
    
    
    eval(strcat("h",num2str(t),"=plot(x,f,'Color',ct);"));
    plot(x,flow,':','Color',ct)
    plot(x,fup,':','Color',ct)

end

xlabel('Growth time (Hours)')
ylabel('P')
legend(Treatments,'Location','northeast');
legend([h1 h2 h3 h4],Treatments,'Location','northeast');  % Only the blue and green lines appear
                                    %   in the legend
ylim([0 1])

grid on

saveas(gcf,fullfile(ResultsFolder,'SurvivalTotal.tif'),'tiffn')  
    
%% PLOT 4

figure (4); % Shape plots
set(gcf, 'Position', get(0, 'ScreenSize')); % Maximize the figure. 
set(gcf,'name','Shape Indexes vs time') 

indexes={'Eccentricity','Solidity','Circularity','# Clusters'};

h1 = plot(rand(1,10));      % Blue line
hold on;
h2 = plot(rand(1,10),'r');  % Red line
h3 = plot(rand(1,10),'g');  % Green line
legend([h1 h3],{'hello','world'});  % Only the blue and green lines appear
                                    %   in the legend
                                    
for index=1:4
     if index==3
        c=-1;        
    else
        c=1;       
    end
    for t=1:4
        subplot(4,4,(index-1)*4+t)
        title(Treatments{t});
        
        Data=reshape(Output{t}(8+index,1:450,:),[450,20]);

        hold on

        for s=1:length(Data(1,:))
            plot(Data(:,s).^c,':')
        end

        plot(nanmean(Data.^c,2))

        xlabel('Frame')
        ylabel(indexes{index})
        xlim([0 450])
        %ylim([0 1])
    end
end

saveas(gcf,fullfile(ResultsFolder,'JOINED_SHAPE_INDEXES.tif'),'tiffn')
%% PLOT 5-8

for fig=5:8
    figure (fig); % Spatial plots
    % Maximize the figure. 
    set(gcf, 'Position', get(0, 'ScreenSize')); 
    set(gcf,'name','Number of clusters over time') 
    
    if fig==7
        c=-1;        
    else
        c=1;       
    end
            
    for t=1:4
        subplot(2,2,t)
        title(Treatments{t});
        Data=reshape(Output{t}(fig+4,1:450,:),[450,20]);

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


function [legend_h,object_h,plot_h,text_strings] = legappend(newStrings,varargin)

h =  findobj(gcf,'Type','axes','Tag','legend');

prop.boxon = get(h,'visible');
prop.loc = get(h,'location'); 
prop.color = get(h,'color'); 
prop.orient = get(h,'Orientation'); 



allDatah = flipud(get(gca,'children')); 
str = get(h,'String'); 

if exist('varargin','var') 
    newStrings = [newStrings,varargin];
end
deleteEntries = sum(cellfun('isempty',newStrings));
if isempty(newStrings) 
    deleteEntries = 1; 
end

if ~deleteEntries
    if iscell(newStrings)
        for k = 1:length(newStrings) 
            str{end+1}=newStrings{k}; 
        end
    end 
    if ~iscell(newStrings)
        str{end+1}=newStrings; 
    end


    [legend_h,object_h,plot_h,text_strings] = legend(h,allDatah,str);

    if strcmpi({prop.boxon},'off')
        legend boxoff
    end

    set(legend_h,'location',prop.loc,'color',prop.color,'Orientation',prop.orient)


end

if deleteEntries
    set(h,'String',str(1:end-nargin))
    [legend_h,object_h,plot_h,text_strings] = legend;
end


if nargout==0
    clear legend_h object_h plot_h text_strings
end


end
