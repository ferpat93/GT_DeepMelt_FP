function [edges]=updateEdges(edges,qEdges,gamma)

mD=mean(edges(edges(:,2)>0,2)); % Mean value of D to scale change f(Q)
qm=mean(qEdges(qEdges>0));
edges(:,3)=max((0.95.*edges(:,3)+((abs(qEdges./qm)).^gamma./(1+abs(qEdges./qm).^gamma)).*mD)./L,0); % Computes new Dij values

end
