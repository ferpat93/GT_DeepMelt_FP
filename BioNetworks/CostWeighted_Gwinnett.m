% Network Cost Gwinnett

% Weighted congestion

path='/Users/lfp3/Dropbox (GaTech)/GT/Spring-19/RootsPaper/ATL_map_files/Gwinett_WT_Dens_Centroids_N5.csv';
%path  = 'C:\Users\lfp3\Dropbox (GaTech)\GT\Spring-19\RootsPaper\ATL_map_files\Gwinett_WT_Dens_Centroids_N5.csv';
Points=dlmread(path,',',1,0);
%w0 = [0 479723 sum(Points(2:end,3)) 3639843-sum(Points(2:end,3))];
%w0 = [0 sum(Points(2:end,3)) 3639843-sum(Points(2:end,3))];
w0=linspace(0,3639843-sum(Points(2:end,3)),20)';

nw=numel(w0);
weights_G = Points(2:end,3);

points=Points(:,1:2); % a has to have just coordinates

% Steiner Network:
gg=load('Gwinnett_SteinerGraph_N5.mat');
GS=gg.GS;
% Leaf Network
gg=load('Gwinnett_LeafGraph_N5_2.mat');
GL=gg.GL;

Gs={GS,GL};
Labels={'Steiner','Leaf'};


EC = cell(2,nw);
Cost = zeros(2,nw);
% Row 1 - ST / Row 2 - LV
%Cols: Population, Area, Length

GIS=[231496 199.38890609 47.185;254594 207.81981171 50.052];
%minC=10000;
%maxC=0;
minC=0.475; % 0.5 in reality
maxC=2.011;

% [ha, pos] = tight_subplot(2,nw,[.01 .01],[.01 .05],[.05 .1]);

for wi=1:nw % Iterate through weights
    
    weights = [w0(wi);weights_G];
    
    
    for n=1:2 % iterate through networks
        G=Gs{n};

        AP=find(or(G.Nodes.Type==1,G.Nodes.Type==3));
        pG = [G.Nodes.X G.Nodes.Y];

        edges_used=[];

        for s=1:length(AP)
            for t=s+1:length(AP)
               [P,d,eu] = shortestpath(G,AP(s),AP(t)); 
               ew = [eu' repmat((weights(s)+weights(t))/2,numel(eu),1)];
               edges_used=[edges_used; ew];
            end
        end

        edge_count = zeros(size(G.Edges,1),1);

        for e=1:size(G.Edges,1)
            edge_count(e) = sum(edges_used(edges_used(:,1)==e,2));
        end

        kk = edge_count./sum(weights);
        minC = min(min(kk),minC);
        maxC = max(max(kk),maxC);
        Cost(n,wi) = sum(G.Edges.Weight.*kk)/1000;
        EC{n,wi} = kk;
        
        %figure(251)
        %subplot(2,nw,wi+nw*(n-1))
%         axes(ha(wi+nw*(n-1)));
%         
%         pe = plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'NodeLabel',{},'EdgeAlpha',1,'LineWidth',kk);
%         pe.EdgeCData = kk;
%         pe.MarkerSize=1;
%         colormap(flip(bone,1));
%         set(gca,'XTick',[]);
%         set(gca,'YTick',[]);
%         caxis manual
%         caxis([minC maxC]);
        %title('Steiner Tree')
    end
end

% Cost Plot

figure(123)
subplot(1,2,1)
plot(w0./pop_gwin,Cost(1,:),'--k')
hold on
plot(w0./pop_gwin,Cost(2,:),'r')
sum(Points(2:end,3))
title('Total Network Cost')
xlabel('Normalized W_o')
ylabel('Cost')
xticks([0 1 2 3])

pop_gwin = sum(Points(2:end,3));
GLL = sum(GL.Edges.Weight)/1000;
GSL = sum(GS.Edges.Weight)/1000;


subplot(1,2,2)
plot(w0./pop_gwin,Cost(1,:)./GSL,'--k')
hold on
plot(w0./pop_gwin,Cost(2,:)./GLL,'r')
sum(Points(2:end,3))
title('Network Cost per unit length')
xlabel('Normalized W_o')
ylabel('Cost per km')
xticks([0 1 2 3])


function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering
if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end
if numel(gap)==1; 
    gap = [gap gap];
end
if numel(marg_w)==1; 
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1; 
    marg_h = [marg_h marg_h];
end
axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
py = 1-marg_h(2)-axh; 
% ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);
    
    for ix = 1:Nw
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
if nargout > 1
    pos = get(ha,'Position');
end
ha = ha(:);
end
