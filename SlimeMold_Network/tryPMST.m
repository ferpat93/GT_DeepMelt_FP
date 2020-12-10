
% Test prunned Minimum spanning tree

s = [1 1 1 2 2 3 3 4 5 5 6 7];
t = [2 4 8 3 7 4 6 5 6 8 7 8];
weights = [10 10 1 10 1 10 1 1 12 12 12 12];
names = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
G = graph(s,t,weights,names);
plot(G,'EdgeLabel',G.Edges.Weight)

T = minspantree(G);


AP = boundary(G.Nodes.X,G.Nodes.Y,0.5);
AP = find(G.Nodes.Type>0);

counts = zeros(1,T.numnodes);

for s=1:length(AP)
    for t=s+1:length(AP)
       P = shortestpath(T,AP(s),AP(t)); 
       counts = counts + histc(P,1:T.numnodes);
    end
end

H = rmnode(T,find(counts==0));

figure;
plot(H)

