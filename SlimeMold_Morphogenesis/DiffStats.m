% Statistical Analysis
Entities={'SlimeMold','Agar','Residuum'};

StructSpotData=load('C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS\OutputData_SpotAnalysis.mat');
SpotData=StructSpotData.Output;

Treatments_Spot={'Glucose_100mM_NaCl_200mM','Glucose_200mM_NaCl_200mM','Glucose_100mM','Glucose_200mM'};
SpotResultsFolder='C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS\';

StructHomoData=load('C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS\RESULTS\OutputData_HomoAnalysis.mat');
HomoData=StructHomoData.Output;

Treatments_Homo={'Glucose_100mM','Glucose_200mM','Control','NaCl_100mM'};
HomoResultsFolder='C:\Users\lfp3\Dropbox (GaTech)\Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS\RESULTS\';


%% 1.Exploration start

ll=0;
ul=1;
nd=200;

xr=linspace(ll,ul,nd);
sp=(ul-ll)/(nd-1);

TT=0.03; % Threshold defined as starting point of exploration
colExp=3; % Column with residuum data

pdT=cell(4,1);

for T=1:4
    Data=reshape(SpotData{T}(colExp,1:450,:),[450,20]); % Get data for the treatment and parameter
    val=zeros(20,1);
    for rep=1:20
        Data(:,rep)
        [fr, ~]= find(Data(:,rep)>TT,1);
        val(rep)=fr*(5/60); % Values in hours
    end
    
    pd = fitdist(val,'Kernel','BandWidth',5); % Get distribution of data
    pdT{T} = pdf(pd,xr);
    
    figure(T)
    histfit(val,5,'kernel')
    hold on
    plot(xr,pdT{T},'k-','LineWidth',2)
    
end

prob=zeros(4,4);

for i=1:4
    for j=1:4
        yi=pdT{i};
        yj=pdT{j};
        diff=max(yi-yj,0);
        prob(i,j)=100*trapz(sp,diff);
    end
end

disp(prob)

%% 2. Bubble Break

ll=0;
ul=15;
nd=200;

xr=linspace(ll,ul,nd);
sp=(ul-ll)/(nd-1);

TT=0.9; % Threshold defined as Solidity for breaking the bubble
colExp=10; % Column with solidity data

pdT=cell(4,1);

for T=1:4
    Data=reshape(SpotData{T}(colExp,1:450,:),[450,20]); % Get data for the treatment and parameter
    val=zeros(20,1);
    for rep=1:20     
        [fr, ~]= find(Data(:,rep)<TT,1);
        val(rep)=fr*(5/60); % Values in hours
    end
    
    pd = fitdist(val,'Kernel','BandWidth',3); % Get distribution of data
    pdT{T} = pdf(pd,xr);
    
    figure(T+5)
    histfit(val,5,'kernel')
    hold on
    plot(xr,pdT{T},'k-','LineWidth',2)
    
end

prob=zeros(4,4);

for i=1:4
    for j=1:4
        yi=pdT{i};
        yj=pdT{j};
        diff=max(yi-yj,0);
        prob(i,j)=100*trapz(sp,diff);
    end
end

disp(prob)

%% 3. Food Reach

ll=5;
ul=40;
nd=200;

xr=linspace(ll,ul,nd);
sp=(ul-ll)/(nd-1);

TT=0; % Threshold defined as starting point of exploration
colExp=7; % Column with residuum data

pdT=cell(4,1);
    
for T=1:4
    Data=reshape(SpotData{T}(colExp,1:450,:),[450,20]); % Get data for the treatment and parameter
    val=zeros(20,1);
    for rep=1:20     
        [fr, ~]= find(Data(:,rep)==TT,1);
        val(rep)=fr*(5/60); % Values in hours
    end
    
    pd = fitdist(val,'Kernel','BandWidth',3); % Get distribution of data
    pdT{T} = pdf(pd,xr);
    
    figure(T+10)
    histfit(val,5,'kernel')
    hold on
    plot(xr,pdT{T},'k-','LineWidth',2)
    
end

prob=zeros(4,4);

for i=1:4
    for j=1:4
        yi=pdT{i};
        yj=pdT{j};
        diff=max(yi-yj,0);
        prob(i,j)=100*trapz(sp,diff);
    end
end

disp(prob)
