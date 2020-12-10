
cases={'A','B','C','D'};
[pathname,FolderSource] = uigetfile('*.txt','Chose Any File inside folder with domains');  % Opens Dialog to select txt that contains Nodes
DomainName={'10-10-1'};
fn=zeros(16,1);

for D=1:length(DomainName)
    FolderRoute=char(strcat(FolderSource,char(DomainName(D)),'/'));
    Case=0;
    for CaseNode=1:1
        for CaseEdge=1:1
            Case=Case+1;
            stringNode=char(strcat(char(cases(CaseNode)),'-Nodes-',char(DomainName(D)),'.txt'));
            stringEdge=char(strcat(char(cases(CaseEdge)),'-Edges-',char(DomainName(D)),'.txt'));
            CaseString=char(strcat('Case',num2str(Case),'-',DomainName(D),'.txt'));
            [fn(Case)]=FlowRoutineFunction(FolderRoute,stringEdge,stringNode,CaseString);
        end
    end
end

%ind=ones(length(HA(:,1)),1);

%for i=1:5
%figure (i)
%histogram(HA(:,i));

%%hold on
%scatter(i.*ind,HA(:,i))
%hold off
%end

%for i=1:length(HA(:,1))
%for i=1:1
%    figure (7)
%    hold on
%    scatter(celledges{i}(:,1),celledges{i}(:,2))
%    hold off
%    figure (8)
%    histogram(celledges{i}(:,1),'Normalization','pdf')
%    figure (9)
%    histogram(celledges{i}(:,1),'Normalization','cdf')
%    figure (10)
%    histogram(celledges{i}(:,2),'Normalization','pdf')
%    figure (11)
%    histogram(celledges{i}(:,2),'Normalization','cdf')
%end
