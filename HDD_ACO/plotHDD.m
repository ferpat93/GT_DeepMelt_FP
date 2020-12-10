function [CoordHDD] = plotHDD(Alignment,HDD_Entry,color)

    hold on
    
    LW=1;
    
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
                plot(CoordHDD([p-1 p],1),CoordHDD([p-1 p],2),'Color',color, 'LineWidth', LW)
                
            else % Vertical or Compound
                
                if Curves(i,5)==0 % Vertical
                    
                    X = Curves(i,3)*(sind(Curves(i,2))-sind(Curves(i,1))) * sign(Curves(i,2)-Curves(i,1));
                    Y = Curves(i,3)*(cosd(Curves(i,1))-cosd(Curves(i,2))) * sign(Curves(i,2)-Curves(i,1));
                    CoordHDD(p,:) = CP + [X Y];
                    
                    plotArc(CoordHDD(p-1,:),Curves(i,:),color)
                    
                else % Compound
                    % Not considered 
                end
                
            end
            
        else % Even -> Straight
            
            CoordHDD(p,:) = CP + Lines(i,1).*[1 tand(Lines(i,2))];
            %figure(10)
            %hold on
            plot(CoordHDD([p-1 p],1),CoordHDD([p-1 p],2),'Color',color, 'LineWidth', LW)
        end

    end

    %figure(10)
    %hold on
    %scatter(CoordHDD(:,1),CoordHDD(:,2))
end



function [] = plotArc(PI,C,color)
    LW=1;
    % Curve: [alpha_v, beta_v, Rv, Rh, horizontal angle]
    t = linspace(C(1),C(2),300);
    r = C(3);
    
    Points = zeros(300,2);
    Points(1,:)= PI;
    
    for i=2:300
        Points(i,:) = Points(i-1,:) + sign(t(i)-t(i-1))*r.*[sind(t(i))-sind(t(i-1)) cosd(t(i-1))-cosd(t(i))];
    end
    
    plot(Points(:,1), Points(:,2),'Color',color, 'LineWidth', LW); 
   
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