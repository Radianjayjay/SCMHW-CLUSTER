%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Statistic MHW of Time, Square, and Mean_Intensity Based on r360x180_60S60N
%%%
%%% In this program, you need to provide the clustering results of the
%%% previous DBSCAN step, the MHWs three-dimensional data and its
%%% corresponding three-dimensional mean intensity data. Finally, you will
%%% get the time, area, cumulative area, intensity, and centroid position
%%% of this data set SCMHW.
%%%
%%% 2024/01/16
%%% @author: Radian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
clc;clear;

% Read the data and store the average intensity in the recognition result [x,y,z,num_idx,mean_int]
disp('Loading Data');

% Enter the file name
filename='idx_DBSCAN_19822022_r360x180_60S60N_Eps5_MinPts150.mat';
filepath=pwd;
mean_int=struct2array(load([filepath '/In/MHW_mean_inXYZ_19822022_r360x180_60S60N.mat']));% % Read mean_int original file

data_idx=load([filepath '/In/' filename]);% Read DBSCAN recognition result file
variable_names = fieldnames(data_idx);
idx=data_idx.(variable_names{1});
location=data_idx.(variable_names{2});

location_idx=[location,idx,mean_int];% [x,y,z,num_idx,mean_int]
clear filename data_idx variable_names idx location mean_int

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store the MHW classification identified by DBSCAN into the structure each_MHW_DBSCAN
% Loop through all events and store them in the structure
disp('Saving each MHW Of DBSCAN in Struct');
each_MHW_DBSCAN=struct();% Creating a structure
MHW_max_num=max(location_idx(:,4));% Count the total number of MHWs identified by DBSCAN
MHW_max_num=200;

for i=1:MHW_max_num

    if mod(i, 10) == 0
        disp(['Separating Rate: ' num2str(i) '/' num2str(MHW_max_num)]);
    end
    
    % Find the row information of all points of the i-th type MHW
    [num_x]=find(location_idx(:,4)==i);
    size_num_x=size(num_x,1);
    
    % Loop to store the [x, y, z, num, mean_int] information of the jth point of the i-th MHW
    for i1=1:size_num_x
        nametocode_i1=['each_MHW_DBSCAN.MHW_DBSCAN_Result.MHW_DBSCAN_' num2str(i) '(' num2str(i1) ',:)=location_idx(num_x(' num2str(i1) '),:);'];
        eval(nametocode_i1);
    end

end

clear location_idx num_x size_num_x i i1 nametocode_i1

% Statistics of the time, area, cumulative area, intensity, and centroid position of each MHW event, Time/Duration/Square/Square_Cum/Mean_Intensity/Centroid
disp('Saving Time/Duration/Square/Square_Cum/Mean_Intensity/Centroid in Struct');
for j=1:MHW_max_num

    % Output progress bar
    if mod(j, 10) == 0
        disp(['Saving T/D/S/Sc/I/C Rate: ' num2str(j) '/' num2str(MHW_max_num)]);
    end
    
    % Calculate the start and end time of the current MHW event and store them in the structure each_MHW_DBSCAN.Time. The kth row represents the time length of the i-th MHW event.
    nametocode_j=['MHW_day_Min=min(each_MHW_DBSCAN.MHW_DBSCAN_Result.MHW_DBSCAN_' num2str(j) '(:,3));MHW_day_Max=max(each_MHW_DBSCAN.MHW_DBSCAN_Result.MHW_DBSCAN_' num2str(j) '(:,3));'];
    eval(nametocode_j);
    each_MHW_DBSCAN.Time(j,1)=(MHW_day_Max-MHW_day_Min+1);
    each_MHW_DBSCAN.Day_Start_End(j,1)=MHW_day_Min;
    each_MHW_DBSCAN.Day_Start_End(j,2)=MHW_day_Max;

    % Calculate and store the area Square and mean intensity Mean_Intensity of the i-th MHW event every day
    for j1=MHW_day_Min:MHW_day_Max
        
        nametocode_j11=['[day_of_square,day_of_square_y]=find(each_MHW_DBSCAN.MHW_DBSCAN_Result.MHW_DBSCAN_' num2str(j) '(:,3)==j1);'];
        eval(nametocode_j11);
        clear day_of_square_y
        sizeday_of_square=size(day_of_square);
        nametocode_j12=['each_MHW_DBSCAN.Square.MHW_DBSCAN_Square_' num2str(j) '(' num2str(j1-MHW_day_Min+1) ',1)=sizeday_of_square(1,1);'];
        eval(nametocode_j12);
        
        % Find the cumulative area Square_Cum
        nametocode_j13=['each_MHW_DBSCAN.Square_Cum.MHW_DBSCAN_Square_Cum_' num2str(j) '(' num2str(j1-MHW_day_Min+1) ',1)=sum(each_MHW_DBSCAN.Square.MHW_DBSCAN_Square_' num2str(j) '(:,1));'];
        eval(nametocode_j13);
        
        % Initialize and find Mean_Intensity
        clear Mean_Intensity_temp
        Mean_Intensity_temp=[];
        
        for j2=1:sizeday_of_square(1,1)
            
            nametocode_j21=['Mean_Intensity_temp(' num2str(j2) ',1)=each_MHW_DBSCAN.MHW_DBSCAN_Result.MHW_DBSCAN_' num2str(j) '(day_of_square(' num2str(j2) '),5);'];
            eval(nametocode_j21);
            
        end
        
        today_mean_int=mean(Mean_Intensity_temp);
        nametocode_j14=['each_MHW_DBSCAN.Mean_Intensity.MHW_DBSCAN_Mean_Intensity_' num2str(j) '(' num2str(j1-MHW_day_Min+1) ',1)=today_mean_int;'];
        eval(nametocode_j14);

        % Find the center of mass at each moment and save it
        Centroid_Tmep=[];
        nametocode_j15=['Centroid_Tmep(:,1:2)=each_MHW_DBSCAN.MHW_DBSCAN_Result.MHW_DBSCAN_' num2str(j) '(day_of_square,1:2);' ...
            'Centroid_Tmep(:,3)=each_MHW_DBSCAN.MHW_DBSCAN_Result.MHW_DBSCAN_' num2str(j) '(day_of_square,5);'];
        eval(nametocode_j15);
        MHW_Mean_Intensity_Mass = sum(Centroid_Tmep(:,3));        
        Centroid_Now_x = sum(Centroid_Tmep(:,1) .* Centroid_Tmep(:,3)) / MHW_Mean_Intensity_Mass;
        Centroid_Now_y = sum(Centroid_Tmep(:,2) .* Centroid_Tmep(:,3)) / MHW_Mean_Intensity_Mass;        
        Centroid_Now = [Centroid_Now_x, Centroid_Now_y];
        nametocode_j16=['each_MHW_DBSCAN.Centroid.MHW_DBSCAN_Centroid_' num2str(j) '(' num2str(j1-MHW_day_Min+1) ',:)=Centroid_Now;'];
        eval(nametocode_j16);

        clear Centroid_Tmep
        clear MHW_Mean_Intensity_Mass Centroid_Now Centroid_Now_x Centroid_Now_y
    end
    
    clear MHW_day_Min MHW_day_Max

end

clear i i1 j j1 j2 
clear nametocode_j nametocode_j11 nametocode_j12 nametocode_j13 nametocode_j14 nametocode_j15 nametocode_j16 nametocode_j21
clear day_of_square sizeday_of_square Mean_Intensity_temp today_mean_int

% Calculate the moving direction Direction, moving speed Speed, and moving distance Distance
% The calculation of moving direction Direction and moving speed Speed ​​is based on the centroid coordinates of each MHW cluster at each moment
for i=1:MHW_max_num

    % Speed
    nametocode_i1=['Temp_Speed_Direction=each_MHW_DBSCAN.Centroid.MHW_DBSCAN_Centroid_' num2str(i) ';'];
    eval(nametocode_i1);
    
    Speed_Displacements = diff(Temp_Speed_Direction, 1, 1);
    Speed_Distance=sqrt(sum(Speed_Displacements.^2, 2));
    Speed_Total_Displacement = sum(Speed_Distance);
    Speed_Total_Time = size(Temp_Speed_Direction,1);
    MHW_Speed_Now = Speed_Total_Displacement / Speed_Total_Time;
    each_MHW_DBSCAN.Speed(i,1)=MHW_Speed_Now;
    
    nametocode_i1=['each_MHW_DBSCAN.Distance.MHW_DBSCAN_Distance_' num2str(i) '=Speed_Distance;'];
    eval(nametocode_i1);
    each_MHW_DBSCAN.Total_Distance(i,1)=Speed_Total_Displacement;
    
    % Maximum range of longitude and latitude in centroid case
    each_MHW_DBSCAN.XY_Range.X_Max_Centroid(i,1)=max(Temp_Speed_Direction(:,1));
    each_MHW_DBSCAN.XY_Range.X_Min_Centroid(i,1)=min(Temp_Speed_Direction(:,1));
    each_MHW_DBSCAN.XY_Range.Y_Max_Centroid(i,1)=max(Temp_Speed_Direction(:,2));
    each_MHW_DBSCAN.XY_Range.Y_Min_Centroid(i,1)=min(Temp_Speed_Direction(:,2));
    
    % The maximum range of longitude and latitude for all points
    nametocode_i2=['Temp_XYZ_Result=each_MHW_DBSCAN.MHW_DBSCAN_Result.MHW_DBSCAN_' num2str(i) ';'];
    eval(nametocode_i2);
    each_MHW_DBSCAN.XY_Range.X_Max_All(i,1)=max(Temp_XYZ_Result(:,1));
    each_MHW_DBSCAN.XY_Range.X_Min_All(i,1)=min(Temp_XYZ_Result(:,1));
    each_MHW_DBSCAN.XY_Range.Y_Max_All(i,1)=max(Temp_XYZ_Result(:,2));
    each_MHW_DBSCAN.XY_Range.Y_Min_All(i,1)=min(Temp_XYZ_Result(:,2));
    
    clear Speed_Displacements Speed_Distance Speed_Total_Displacement Speed_Total_Time MHW_Speed_Now nametocode_i1
    clear nametocode_i1 nametocode_i2 Temp_XYZ_Result

    % Direction
    Halfway_Point = floor(size(Temp_Speed_Direction, 1) / 2);
    Temp_Speed_Direction_Satrt = Temp_Speed_Direction(1:Halfway_Point, :);
    Temp_Speed_Direction_End = Temp_Speed_Direction(Halfway_Point+1:end, :);
    Average_Satrt = mean(Temp_Speed_Direction_Satrt, 1);
    Average_End = mean(Temp_Speed_Direction_End, 1);
    Direction_Angle = atan2(Average_End(2) - Average_Satrt(2), Average_End(1) - Average_Satrt(1));
    Direction_Angle_Degrees = rad2deg(Direction_Angle);
    each_MHW_DBSCAN.Direction(i,1)=Direction_Angle_Degrees;
    
    % Square_Max,Square_Cum_Max
    nametocode_i3=['Temp_Square_Max=each_MHW_DBSCAN.Square.MHW_DBSCAN_Square_' num2str(i) ';'];
    eval(nametocode_i3);
    nametocode_i4=['Temp_Square_Cum_Max=each_MHW_DBSCAN.Square_Cum.MHW_DBSCAN_Square_Cum_' num2str(i) ';'];
    eval(nametocode_i4);
    each_MHW_DBSCAN.Square_Max(i,1)=max(Temp_Square_Max);
    each_MHW_DBSCAN.Square_Cum_Max(i,1)=max(Temp_Square_Cum_Max);

    % Mean_Intensity_Total
    nametocode_i5=['Temp_Mean_Intensity_Total=each_MHW_DBSCAN.Mean_Intensity.MHW_DBSCAN_Mean_Intensity_' num2str(i) ';'];
    eval(nametocode_i5);
    each_MHW_DBSCAN.Mean_Intensity_Total(i,1)=mean(Temp_Mean_Intensity_Total);

    clear Halfway_Point Temp_Speed_Direction Temp_Speed_Direction_Satrt Temp_Speed_Direction_End Average_Satrt Average_End Direction_Angle Direction_Angle_Degrees
    clear nametocode_i3 nametocode_i4 nametocode_i5 Temp_Square_Max Temp_Square_Cum_Max Temp_Mean_Intensity_Total

end

clear i MHW_max_num

% Store the final result structure: each_MHW_DBSCAN
disp('Saving struct: each_MHW_DBSCAN');

% Store the output file in the /Out folder according to the path information
save([filepath '/Out/each_MHW_DBSCAN_Eps5_MinPts150.mat'],'each_MHW_DBSCAN','-v7.3');