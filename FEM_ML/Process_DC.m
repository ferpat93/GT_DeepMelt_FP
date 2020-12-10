% Process Disp controlled data

close all
clc

root_folder = pwd;
root_folder = 'C:\Users\lfp3\Documents\Local_AutoAbaqus\';

folder = fullfile(root_folder,'Disp_OutFiles');
Processed_folder = fullfile(root_folder,'Processed_Output');

sims = GetFiles(folder,'.out');

InputPar = readtable(fullfile(pwd,'Disp_InputParameters.txt'),'FileType','text');
InputPar.sv = (9.81/1000).*InputPar.d.*InputPar.Depth; %kPa
InputPar.sh = InputPar.sv.*(1-sind(InputPar.phi)); %kPa

rc = 0.5; % radius of cavity
c = 5 ; % Cohesion kPa


InputPar.Delft_inf = GetDelft(InputPar,rc,inf,c); % (T,r,rmax,c)
InputPar.Delft_sug = GetDelft(InputPar,rc,0.67.*InputPar.Depth,c); % (T,r,rmax,c)
InputPar.Exists = false([size(InputPar,1) 1]);
InputPar.yield_angle = nan([size(InputPar,1) 1]);
InputPar.yield_stress = inf([size(InputPar,1) 1]);
InputPar.yield_strain = inf([size(InputPar,1) 1]);

%Wblock = GetBlockW(InputPar,rc); % Weight of block
%Yield_Elastic = GetElastic(InputPar,c,angles);
%{
figure(12222); hold on
for s=1:numel(sims)
    plot(angles,Yield_Elastic(s,:));
end
%}

for si=1:numel(sims)
    s = find(strcmp(InputPar.Folder,sims{si}(1:end-4)));
    InputPar.Exists(s) = true;
    
    sv = InputPar.sv(s); % Vertical Stress
    sh = InputPar.sh(s); % Horizontal Stress

    Data = readmatrix(fullfile(folder,sims{si}),'FileType','text');
    
    Ps = unique(Data(:,1));
    np = size(Data,1)/numel(Ps);
    keys = 1:np:size(Data,1);   
    angles = linspace(-90,90,np);
      
    [~,order] = sortrows(Data(1:np,:),3);
    Data3 = zeros(np,numel(keys),size(Data,2));
    for i=1:numel(keys)
        Data3(:,i,:) = Data(order+(keys(i)-1),:);
    end
    clear Data

    Data3([1 np],:,6) = Data3([1 np],:,6).*2; Data3([1 np],:,7) = Data3([1 np],:,7).*2;
    Data3(:,:,8) = (Data3(:,:,6).^2 + Data3(:,:,7).^2).^(1/2); % Add resulting force
    Data3(:,:,9) = Data3(:,:,8)./(1000*rc*deg2rad(180/(np-1))); % Add eq. stress (kPa)
    Data3(:,:,10) = (Data3(:,:,4).^2 + Data3(:,:,5).^2).^(1/2); % Add displacement
      
    e = Data3(1,:,10)/rc; % radial strain
    Stress = Data3(:,:,9); % Tangential stress
    Out = [0 angles ; e' Stress'];
    writematrix(Out,fullfile(Processed_folder,sims{si}),'FileType','text');
    
%     [X,Y] = meshgrid(e,angles);
%     figure(1000+s); surf(X,Y,Stress)
  
    % Find yield points
    peaks = nan(np,4);
    for a=1:np
        [p,l]=findpeaks(Stress(a,:),'NPeaks',1);
        if ~isempty(p); peaks(a,:) = [p angles(a) e(l) a]; % radial strain - peak stress
        end
    end
    peaks(isnan(peaks(:,1)),:) = [];
    peaks = sortrows(peaks,3);
    
    if ~isempty(peaks)    
        InputPar.yield_angle(s) = peaks(1,2);
        InputPar.yield_stress(s) = peaks(1,1);
        InputPar.yield_strain(s) = peaks(1,3);
        %under_delft = find(peaks(:,1)<Delft(s));
        %cols = [0 0 0; 1 0 0];
        %figure(s+3000);
        %scatter(peaks(:,2),peaks(:,3),90.*normalize(peaks(:,1),'range')+10,cols(1+(peaks(:,1)<Delft(s)),:),'filled')
        %to_plot = [round(linspace(1,np,3)) peaks(under_delft(1),4)]; % Angles to plot
    %else        
        
    %    to_plot = round(linspace(1,np,3)); % Angles to plot     
    end
    
%    to_plot = round(linspace(1,np,3));
%     
%     figure(si+1000); hold on
%     for i=1:numel(to_plot)
%         S = Stress(to_plot(i),:);
%         %subplot(1,3,3); hold on
%         plot(e,S)
%     end
%        
%     plot([0 e(end)],[sv sv]); plot([0 e(end)],[sh sh]);% plot([0 e(end)],[InputPar.Delft(s) InputPar.Delft(s)])
%     annotation('textbox',[.9 .5 .1 .2],'String',['H = ' num2str(InputPar.Depth(s))],'EdgeColor','none')
%         
end

writetable(InputPar,'Disp_Output.csv')
T = InputPar;

T = readtable('Disp_Output.csv');
T.sinphi = sind(T.phi);
T.sinpsi = sind(T.psi);

TT = T;
TT.Folder = [];
M = table2array(TT);
N = normalize(TT);

figure; hold on
sz = 20;
scatter(InputPar.Delft_sug./1000,InputPar.yield_stress./1000,sz,1-sind(InputPar.phi),'filled')
plot([0 10],[0 10],'r')
colormap('copper') ; colorbar
xlabel('P_{lim} [MPa] - Delft Equation'); ylabel('Yield Stress [MPa]')

%% Auxiliar Functions 

function [R] = GetSymmetry(Stress)

    A = Stress(2:18,:);
    Sr = sum((A-flipud(Stress(20:end-1,:))).^2);
    St = sum((A-mean(Stress(2:18,:))).^2);
    St2 = (max(A)-min(A)).^2;
    C = 1-Sr./St2;
    
end


function [P] = GetDelft(T,r,rmax,c)
    senn = sind(T.phi); coss = cosd(T.phi); cott = coss./senn;
    
    G = T.E./(2000.*(1+T.v)); % Shear mod. kPa
    Q = (T.sv .* senn + c.*coss)./G ;
    P = (T.sv.*(1+senn)+c.*(coss+cott)).*((r./rmax).^2+(Q)).^(-senn./(1+senn)) - c*cott;

end

function [P] = GetElastic(T,c,angles)
    senn = sind(T.phi); K = 1-senn; Kp = (1+senn)./(1-senn);
    Y = 2*c .* (cosd(T.phi)./K); po = T.sv;
    
    A = ((1+K).*po-Y)./(1+Kp); B = (2.*po.*(1-K))./(1+Kp);
    
    P = zeros(size(T,1),numel(angles));
    
    for i =1:numel(angles)
        P(:,i) = A + B.*cosd(2*angles(i));
    end

end


function [W] = GetBlockW(T,r)
    tann = tand(T.phi./2); 
    sv = T.sv;
    
    W = sv.*(T.Depth.*tann+2*r);
end

function [names]=GetFiles(path,string)
        d = dir(path);
        isub = [d(:).isdir]; %# returns logical vector
        names = {d(~isub).name}';
        %names(or(ismember(names,{'.','..'}),~contains(names,string))) = [];
        names(~contains(names,string)) = [];
end 
