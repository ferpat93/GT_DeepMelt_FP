% [CoordHDD] = plotHDD(Alignment,HDD_Entry,color)

% Profile View

close all


Plan = [0 0; 20 20; 30 30];

folder = '/Users/lfp3/Dropbox (GaTech)/GT/Fall_19/HDD-ACO/';
Alignment = readmatrix(fullfile(folder,'AlignmentFig1.csv')); % Assumes Actual design was there 
%Profile = readmatrix(fullfile(folder,'Profile.csv'));
%Profile(end,:)=[];

figure(10)
hold on
subplot(1,2,2)
title('Profile View')
%plot(Profile(:,1),Profile(:,2),'k-')
axis tight

Entry_station = 0;
HDD_Entry = [Entry_station 0];
PV=plotHDD2(Alignment,HDD_Entry);

%ylim([-110 20])
%xlim([0 2050])
xlabel('Station')
ylabel('Elevation')
axis tight


ax = gca;
ax.XRuler.Exponent = 0;



% Plan View

Plan_Angles = [22 22 22 22 38 38 38];
Colors = [0 38 154;0 38 154;150 9 22;150 9 22;150 9 22;0 38 154;0 38 154]./255;
Lines = {'-','-.'};
Coords = zeros(size(Alignment,1)+1,2);

figure(10)
subplot(1,2,1)
hold on
xlabel('X')
ylabel('Y')
title('Plan View')

for s=2:size(Alignment,1)+1
    
    
    if abs(s-5)>0   
        Coords(s,:) = Coords(s-1,:) + Alignment(s-1,1).*[cosd(Plan_Angles(s-1)) sind(Plan_Angles(s-1))];
        
        plot(Coords(s-1:s,1),Coords(s-1:s,2),'Color',Colors(s-1,:), 'LineWidth', 2,'LineStyle',Lines{mod(s,2)+1})
        
    else 
        % Curve: [alpha_v, beta_v, Rv, Rh, horizontal angle]
        X = 1000*(sind(38)-sind(22));
        Y = 1000*(cosd(22)-cosd(38));          
        Coords(s,:) = Coords(s-1,:) + [X Y];
        
        plotArc(Coords(s-1,:),[22 38 1000 0 0],Colors(s-1,:))
    
    end
    
end


figure(10)
subplot(1,2,2)
hold on
scatter(PV(:,1),PV(:,2),20,'k','filled')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

figure(10)
subplot(1,2,1)
hold on
scatter(Coords(:,1),Coords(:,2),20,'k','filled')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
xlim([-50 1100])


function [] = plotArc(PI,C,color)
    LW=2;
    % Curve: [alpha_v, beta_v, Rv, Rh, horizontal angle]
    t = linspace(C(1),C(2),300);
    r = C(3);
    
    Points = zeros(300,2);
    Points(1,:)= PI;
    
    for i=2:300
        Points(i,:) = Points(i-1,:) + sign(t(i)-t(i-1))*r.*[sind(t(i))-sind(t(i-1)) cosd(t(i-1))-cosd(t(i))];
    end
    
    plot(Points(:,1), Points(:,2),'Color',color, 'LineWidth', LW,'LineStyle','-.'); 
   
end
function [Lines, Curves] = TranslateAlignment(Alignment)
    
    nSegments = size(Alignment,1);
    
    Lines = zeros(floor(nSegments/2)+1,2);
    Curves = zeros(floor(nSegments/2),5);
    
    for s = 1:nSegments
        
        i = (s+mod(s,2))/2;
        
        if mod(s,2)==1 % Odd -> Straight
            Lines(i,:)=Alignment(s,1:2); % [Plan Length Angle]
        else % Curve
            Curves(i,:)=[Alignment(s-1,2) Alignment(s+1,2) Alignment(s,4) Alignment(s,6) Alignment(s,5)]; % [alpha_v, beta_v, Rv, Rh, horizontal angle]
        end  
        
    end

end
function [CoordHDD] = plotHDD2(Alignment,HDD_Entry)

    hold on
    
    LW=2;
    Colors = [0 38 154;0 38 154;150 9 22;150 9 22;150 9 22;0 38 154;0 38 154]./255;
    
    [Lines, Curves] = TranslateAlignment(Alignment);
    
    % Lines -> 3, 4 or 5 rows.
    % line i : [plan length, angle]

    % Curves -> 2, 3 or 4 rows.
    % row i : [alpha_v, beta_v, Rv, Rh, horizontal angle]

    % HDD Entry: [Station Elevation]

    % GS: Ground Surface [Stations/Elevations]

    nPoints = size(Lines,1) + size(Curves,1) + 1;

    CoordHDD = zeros(nPoints, 2);
    CoordHDD(1,:) = HDD_Entry;
    
    for p = 2:nPoints
        
        CP = CoordHDD(p-1,:);
        i = (p-mod(p,2))/2;
        
        if mod(p,2)==1 % Odd -> Curve
            PL = Curves(i,4)*deg2rad(Curves(i,5));
            if Curves(i,1)==Curves(i,2) % Horizontal or Helix
                
                CoordHDD(p,:) = CP + PL.*[1 tand(Curves(i,2))];
                plot(CoordHDD([p-1 p],1),CoordHDD([p-1 p],2),'Color',Colors(p-1,:), 'LineWidth', LW,'LineStyle','-.')
                
            else % Vertical or Compound
                
                if Curves(i,5)==0 % Vertical
                    
                    X = Curves(i,3)*(sind(Curves(i,2))-sind(Curves(i,1))) * sign(Curves(i,2)-Curves(i,1));
                    Y = Curves(i,3)*(cosd(Curves(i,1))-cosd(Curves(i,2))) * sign(Curves(i,2)-Curves(i,1));
                    CoordHDD(p,:) = CP + [X Y];
                    
                    plotArc(CoordHDD(p-1,:),Curves(i,:),Colors(p-1,:))
                    
                else % Compound
                    % Not considered 
                end
                
            end
            
        else % Even -> Straight
            
            CoordHDD(p,:) = CP + Lines(i,1).*[1 tand(Lines(i,2))];
            %figure(10)
            %hold on
            plot(CoordHDD([p-1 p],1),CoordHDD([p-1 p],2),'Color',Colors(p-1,:), 'LineWidth', LW)
        end

    end

    %figure(10)
    %hold on
    %scatter(CoordHDD(:,1),CoordHDD(:,2))
end
