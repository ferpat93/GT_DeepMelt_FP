% Analyze indexes EP array

load('Cases.mat') 
Cases(2,:)=[];
% [E,v,phi,psi]

%   [Eccentricity,Solidity,MinorAxis,MajorAxis,Area,MaxDist,AngletoMax]
load('stats_array_EP');
load('stats_array_CAV');
load('contourEP');
statsCAV = Cav_Indexes;

nModels=numel(contourEP);

%EP: [Eccentricity,Solidity,MinorAxis,MajorAxis,Area,MaxDist,AngletoMax]
indexNames_EP = {'Eccentricity','Solidity','MinorAxis','MajorAxis','Area','Max. Dist from origin','Angle to farthest plastic point'};
ind_plot = [5,2,6,7];
% CAV: [+Axis -Axis Ecc Solidity Area Alt_+Axis Alt_-Axis]
indexNames_CAV = {'MinorAxis','MajorAxis','Eccentricity','Solidity','Area','Alt+Axis','Alt-Axis'};


% Classify EP'S

[cido,cmeans] = kmeans(statsEP(:,[1 2 5 6 7]),3,'dist','correlation'); % 'correlation' 'sqeuclidean'
%[silh,h] = silhouette(statsEP(:,[1 2 5 6 7]),cid,'sqeuclidean');

cid = 2*ones(nModels,1);
max_id=find(cmeans(:,3)==max(cmeans(:,3)));
min_id=find(cmeans(:,3)==min(cmeans(:,3)));
cid(cido==max_id)=3;
cid(cido==min_id)=1;

A = [Cases(:,1:2) sind(Cases(:,3:4)) cid];
colors = [0 0 1 ; 0 1 0; 1 0 0];
c=colors(cid,:);


%% Relationship between Indexes

% EP: [Eccentricity,Solidity,MinorAxis,MajorAxis,Area,MaxDist,AngletoMax]
vs = [ 2 7; 2 1; 2 6];
for p=1:size(vs,1)
    VSplot(vs(p,1),vs(p,2),statsEP,indexNames_EP,c)
end

% CAV: [+Axis -Axis Ecc Solidity Area Alt_+Axis Alt_-Axis]
vs = [3 5];
for p=1:size(vs,1)
    VSplot(vs(p,1),vs(p,2),statsCAV,indexNames_CAV,c)
end

% COMPARE EP TO EP INDEXES -> There is no correlation!
scatter(statsCAV(:,5),statsEP(:,6),[],c) % For example
figure
scatter(Cases(:,1)./sind(Cases(:,3)),statsEP(:,5),[],c) % For example
figure
scatter(Cases(:,4),statsEP(:,5),[],c) % For example



%% PLOT EP CONTROUS BY GROUP

for m=1:nModels
    id = cid(m);
    cont = contourEP{m};
    
    figure(1000)
    hold on
    if id==1
        h1=plot(cont(:,1),cont(:,2),'Color',colors(id,:));
    elseif id==2
        h2=plot(cont(:,1),cont(:,2),'Color',colors(id,:));
    else
        h3=plot(cont(:,1),cont(:,2),'Color',colors(id,:));
    end
    
end

legend([h1,h2,h3],'Low','Mid','High')

%% INDEX HISTOGRAMS

% Split histograms EP (Do the same for the mechanical parameters?)
nC = max(cid);

for i=1:numel(ind_plot)
    figure(2000)
    ax2 = subplot(2,2,i);
    hold on
    for cat=1:nC
        histogram(statsEP(cid==cat,ind_plot(i)),'FaceColor',colors(cat,:),'FaceAlpha',0.6) % ,'Normalization','probability'
    end
    if i==1
        set(ax2, 'XScale', 'log')
    end
    title(indexNames_EP{ind_plot(i)})  
end

% Split histograms CAV (Do the same for the mechanical parameters?)
nC = max(cid);
rr=[3,5];
for i=1:2
    figure(3000)
    subplot(1,2,i)
    hold on
    for cat=1:nC
        histogram(statsCAV(cid==cat,rr(i)),'Normalization','pdf','FaceColor',colors(cat,:),'FaceAlpha',0.6)
    end
    title(indexNames_CAV{rr(i)})  
end


%% To this point the classification looks pretty good... now try to involve soil parameters

%{
1. find mf that is missing
2. load mechanical parameters (exclude missing)
3. plot 3d groups (excluding dilation)
4. train classification thing - try with 2
and 3 clusters
%}

% ptsymb = {'bs','r^','md','go','c+'};
% for i = 1:2
%     clust = find(cidx2==i);
%     plot3(meas(clust,1),meas(clust,2),meas(clust,3),ptsymb{i});
%     hold on
% end

%% AUX FUNCTIONS

function []=VSplot(iX,iY,statsEP,indexNames,c)

    figure
    scatter(statsEP(:,iX),statsEP(:,iY),[],c)
    title (strcat(indexNames{iX},' VS. ',indexNames{iY}))
    xlabel(indexNames{iX}) 
    ylabel(indexNames{iY}) 

end