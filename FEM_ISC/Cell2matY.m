YCav=zeros(625,4);
YEP=zeros(625,3);

for i=1:625
    
    YCav(i,1)=((stats{i}.Area/pi())-1)*100;  % Percentage of increase 
    YCav(i,2)=stats{i}.Eccentricity;
    YCav(i,3)=stats{i}.Solidity;
    YCav(i,4)=stats{i}.Alt_MajorAxis/stats{i}.Alt_MinorAxis;
end

for i=1:624
    YEP(i,1)=(statsEP{i}.Area-stats{i}.Area)/pi();  % Percentage of increase 
    YEP(i,2)=statsEP{i}.Eccentricity;
    YEP(i,3)=stats{i}.Solidity;
end


Cfullpath = fullfile('C:\Users\lfp3\Dropbox\GT\Spring-18\IS_Paper\Matlab',filesep,'Y_Responses'); % Gets info of the files inside the given folder with the given extension
save(Cfullpath,'YCav','YEP','-v7.3');