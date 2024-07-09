%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save MHWs in XYZ data Based on r360x180_60S60N
%%%
%%% In this program, you need to provide the Lat, Lon, and Time information
%%% of the data you use in the In folder and store it as a mat file. After
%%% adjusting the parameters, you will get a three-dimensional array containing
%%% the MHW three-dimensional spatial information, and also a three-dimensional
%%% array of the corresponding mean intensity. 
%%%
%%% 2023/12/20
%%% @author: Radian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
clear,clc;

% Read longitude and latitude
disp('Loading: lon_lat_time_num');
load *:\*\02_EveryDayMHWinXYZ\In\lon_lat_time_num_r360x180_60S60N.mat data_time_num;%Windows
% load */*/02_EveryDayMHWinXYZ/In/lon_lat_time_num_r360x180_60S60N.mat data_time_num;%Linux

% Read MHW cell data
disp('Loading: MHW Cell');
load *:\*\02_EveryDayMHWinXYZ\In\MHW_oisst_1982_2022_r360x180_60S60N_detectc.mat MHWc;% Windows
% load /*/02_EveryDayMHWinXYZ/In/MHW_oisst_1982_2022_r360x180_60S60N_detectc.mat MHWc;% Linux

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create NaN array
disp('CreatingNaN');
[NaN_X,NaN_Y]=size(MHWc);
NaN_Z=max(size(data_time_num));
meaninXYZ_19822022_r360x180_60S60N=NaN(NaN_X,NaN_Y,NaN_Z);

% Store the mean_int data in the Cell array into the NaN array point by point
for i=1:NaN_X

    disp(['Saving Lon：' num2str(i) '/' num2str(NaN_X)]);

    for j=1:NaN_Y
        
        % disp(['Saving Lat：' num2str(j) '/' num2str(NaN_Y)]);
        % Temporarily store the current grid data in temp
        temp=MHWc{i,j};
        sizetemp=size(temp);

       for k=1:sizetemp(1,1)

            % disp(['Saving mean：' num2str(k) '/' num2str(sizetemp(1,1))]);
            [timestart,Ystart]=find(temp(k,1)==data_time_num(:,1));
            [timeend,Yend]=find(temp(k,2)==data_time_num(:,1));
            mean_int=temp(k,5);
            meaninXYZ_19822022_r360x180_60S60N(i,j,timestart:timeend)=mean_int;

       end

    end

end

disp('Saving OutDATA');
save('*:\*\02_EveryDayMHWinXYZ\Out\meaninXYZ_19822022_r360x180_60S60N.mat','meaninXYZ_19822022_r360x180_60S60N');% Windows
% save('/*/02_EveryDayMHWinXYZ/Out/meaninXYZ_19822022_r360x180_60S60N.mat','meaninXYZ_19822022_r360x180_60S60N');% Linux

