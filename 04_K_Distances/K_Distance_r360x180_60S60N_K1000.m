%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate K distance Based on r360x180_60S60N
%%%
%%% This program calculates the K distance table of the three-dimensional
%%% data array generated in the previous step to determine the parameters of
%%% the DBSCAN algorithm.
%%%
%%% 2023/12/22
%%% @author: Radian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear,clc;

% Set the calculated K distance range
K_near=1000;
filepath=pwd;
step=[19822022,19821983,19971998,20092010,20152016];
sizestep=size(step);

for i=1:max(sizestep)

    % Read the data file, assuming the data file is data.txt, each line contains x, y, z coordinates
    disp(['Loading: ' num2str(step(i))]);
    X=struct2array(load([filepath '/In/MHW_location_inXYZ_' num2str(step(i)) '_r360x180_60S60N.mat']));
    
    % Calculate the K distance for the specified number of neighbors
    disp(['K-Distances: ' num2str(step(i))]);
    tic
    KD = pdist2(X,X,'euc','Smallest',K_near);
    toc
    
    % Store results
    disp(['Saving: ' num2str(step(i))]);
    save([filepath '/Out/K_Distance_' num2str(step(i)) '_r360x180_60S60N_K' num2str(K_near) '.mat'],'KD','-v7.3');

end