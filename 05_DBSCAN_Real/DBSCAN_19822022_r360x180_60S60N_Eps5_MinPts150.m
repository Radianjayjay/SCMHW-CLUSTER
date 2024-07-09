%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DBSCAN
%%%
%%% This program performs clustering operations on three-dimensional data
%%% by setting appropriate Eps and minPts
%%%
%%% 2023/01/12
%%% @author: Radian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization
clc;clear;
Eps=[5];% Set the distance parameter in DBSCAN
MinPts=[150];% Set the scatter point number threshold in DBSCAN

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loading');
% Read the location information of non-NaN data% Read the location information of non-NaN data
filepath=pwd;
location=struct2array(load([filepath '/In/MHW_location_inXYZ_19822022_r360x180_60S60N.mat']));

sizeEps=size(Eps);
sizeMinPts=size(MinPts);
disp('DBSCAN_ing');
for i=1:sizeEps(2)
    for j=1:sizeMinPts(2)
        fprintf('Eps: %s\nMinPts: %s\n', num2str(Eps(i)), num2str(MinPts(j)));
        tic
        idx_DBSCAN_19822022_r360x180_60S60N = dbscan(location,Eps(i),MinPts(j));
        toc
        disp('Saving');
        save([filepath '/Out/idx_DBSCAN_19822022_r360x180_60S60N_Eps' num2str(Eps(i)) '_MinPts' num2str(MinPts(j)) '.mat'],'idx_DBSCAN_19822022_r360x180_60S60N','location','-v7.3');
    end
end