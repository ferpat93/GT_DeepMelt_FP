% Check Porosity at boundaries

%load('Data1.mat');

[xx,yy,zz] = meshgrid(1:100,1:100,1:84);
[theta,rho,z] = cart2pol(xx(:)-50,yy(:)-50,zz(:));

sizeI = [100 100 84];
T = zeros(prod(sizeI),13,'single');

for i=1:13
    CP = imresize3(Porosity{i,3}-Porosity{i,1},sizeI);
    T(:,i) = abs(CP(:));
end

interval = 5;
range = 0:interval:50;
avg_pc = zeros(numel(range)-1,13,2);

for i=1:numel(range)-1
    S = T(and(rho>=range(i),rho<range(i+1)),:);
    avg_pc(i,:,1) = nanmean(S);
    avg_pc(i,:,2) = prctile(S,90);  
end

xvalues = (range(1:end-1)+range(2:end))./2;

TestVariables = [ [2 4 6 2 4 6 2 4 2 6 6 4 6]' [1 1 1 1 1 1 0 0 0 0 0 0 0]' [0 0 0 1 1 1 0 0 1 0 1 1 2]']; % Length (2,4,6), density (0,1), load(0,1,2)
colorlines = {'r','k','b'};
linetyp = {'-.','-'};


figure;

for l=2:2:6
    
    idx = find(TestVariables(:,1)==l);
    
    subplot(1,3,l/2); hold on
    for i=1:numel(idx)
        plot(xvalues,avg_pc(:,idx(i),2),[colorlines{TestVariables(idx(i),3)+1} linetyp{TestVariables(idx(i),2)+1}])
    end
    
    %ylim([0 11]); 
    xlabel('Distance from center [mm]'); 
    ylabel('Average absolute porosity change [%]'); 
    title(strcat('L =',{' '},num2str(l),'D')); grid on
end

legend({'Dense - \sigma1','Dense - \sigma2','Loose - \sigma1','Loose - \sigma2','Loose - \sigma3'},'Location','southoutside');

