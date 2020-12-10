
[Name,Path] = uigetfile('*.txt','Chose Topology file'); 
topo=dlmread(strcat(Path,Name));     % Stores as 'nodes' the  matrix


[n1,wout1]=hist(topo(:,3));
% Normalize so that the sum of the points = 1.
%normalizedCounts = n1 / sum(n1);
plot(wout1);


%plot count and bincnter 
h1 = 0.05 * randn(1, 10000);
h2 = 0.10 * randn(1, 10000);
h3 = 0.30 * randn(1, 10000);
[counts1, binCenters1] = hist(h1, 500);
[counts2, binCenters2] = hist(h2, 500);
[counts3, binCenters3] = hist(h3, 500);
plot(binCenters1, counts1, 'r-');
hold on;
plot(binCenters2, counts2, 'g-');
plot(binCenters3, counts3, 'b-');
grid on;
% Put up legend.
legend1 = sprintf('mu = %.3f', mean(h1));
legend2 = sprintf('mu = %.3f', mean(h2));
legend3 = sprintf('mu = %.3f', mean(h3));
legend({legend1, legend2, legend3});