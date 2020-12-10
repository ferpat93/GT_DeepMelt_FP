%% Process Cavity

% Assuming every simulation was completed succesfully

Cavity_folder='';

cd(Cavity_folder)

Files=dir('*.csv');

parfor f=1:Length(Files)
    File=importCav(Files(f));
    
end
