% Extract Folders

names=extractfield(subFolders,'name')';
names(1:2)=[];

dates=cell(80,1);

for i=1:length(names)
    a=strsplit(names{i},'-');
    names{i}=a{2};
    a=strsplit(a{2},'_');
    dates{i}=strcat(a{end-3},'_',a{end-2});
end