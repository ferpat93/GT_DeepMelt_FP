global PEDGE PNODE

%% INPUT VARIABLES
NX=40; %Nodes on X direction
NY=50; %Nodes on Y direction
D=1; %Lattice distance (side of equilateral triangle)

NE=3*NX*NY-2*(NX+NY)+1;

PEDGE=zeros(NE,1);
PNODE=zeros(NY,NX);

for i=1:NE
PEDGE(i)=rand();
end

for k=1:NX %Iterates over X
    for f=1:NY %Iterates over Y
    PNODE(f,k)=rand(); %set normally distributed size of pore
    end
end