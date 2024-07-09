%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find MHWs Based on r360x180_60S60N
%%%
%%% In this program, you need to provide the NC data of SST. After adjusting 
%%% the parameters and running the program, you will get a mat file in the 
%%% form of Cell array.
%%%
%%% 2023/12/20
%%% @author: Radian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear,clc;

% Please set the path of m_mhw1.0 toolbox, 
% you can get the toolbox from the following link https://github.com/ZijieZhaoMMHW/m_mhw1.0
addpath('*:\*\m_mhw1.0_202312');%window;
% addpath('/*/m_mhw1.0_202312/');%Linux

% Set the initial file name. Here you can set the path and name of the input SST data NC file.
filename = '*:\*\sst_oisst_1982_2022_r360x180_60S60N.nc';%window;
% filename = '/*/sst_oisst_1982_2022_r360x180_60S60N.nc';%Linux

% Set data length, set climate state range, set MHW identification range
data_start=datenum(1982,1,1);
data_end=datenum(2022,12,31);
clim_start=datenum(1982,1,1);
clim_end=datenum(2022,12,31);
find_MHW_start=datenum(1982,1,1);
find_MHW_end=datenum(2022,12,31);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the file
% Prompt the current operation file name
disp(['Loading:' filename]); 

% Read all variable names in the NC file
info = ncinfo(filename);
varNames = {info.Variables.Name};

% Select the variable name to read
selectedVarName_1 = varNames{1,1};%lon
selectedVarName_2 = varNames{1,2};%lat
selectedVarName_3 = varNames{1,3};%time
selectedVarName_4 = varNames{1,4};%sst

% Use ncread function to read variable data
data_lon = ncread(filename, selectedVarName_1);
data_lat = ncread(filename, selectedVarName_2);
data_time = ncread(filename, selectedVarName_3);
data_sst = ncread(filename, selectedVarName_4);

% Replace the original time stamp
t0 = datenum(1982,1,1);
data_time=t0+data_time-data_time(1,1);

clear filename info varNames selectedVarName_1 selectedVarName_2 selectedVarName_3 selectedVarName_4;
clear data_lon data_lat data_time t0;

% Identify MHW detectc
tic
disp('NOW: Finding MHW'); 
[MHWc,mclimc,m90c,mhw_tsc]=detectc(data_sst,data_start:data_end,clim_start,clim_end,find_MHW_start,find_MHW_end);
toc

% Store data results detectc
disp('Saving:');
save('*:\*\MHW_oisst_1982_2022_r360x180_60S60N_detectc.mat','MHWc','mclimc','m90c','mhw_tsc','-v7.3');%Windows
% save('/*/MHW_oisst_1982_2022_r360x180_60S60N_detectc.mat','MHWc','mclimc','m90c','mhw_tsc','-v7.3');%Linux