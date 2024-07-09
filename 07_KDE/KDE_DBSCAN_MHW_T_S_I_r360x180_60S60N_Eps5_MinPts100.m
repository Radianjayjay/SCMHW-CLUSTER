%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% KDE_DBSCAN_MHW_T_S_I_r360x180_60S60N_Eps5_MinPts150
%%%
%%% In this program, based on the three indicators of time, area and
%%% average intensity counted in the previous step, the KDE algorithm is
%%% used to take the 90% threshold of each indicator as the reference point
%%% to classify all SCMHW into 8 categories.
%%%
%%% 2024/01/18
%%% @author: Radian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The [Time, Square, Mean_Intensity] of the heat wave event clusters were sorted and KDE classification was performed to obtain 8 categories of results
clc;clear;

% Set the data start and end time
Datestart = datenum('1982-01-01');
Dateend = datenum('2022-12-31');

% Set the boundary between extreme events and non-extreme events
KDE_Time_Percentile_Value=0.9;
KDE_Square_Percentile_Value=0.9;
KDE_Mean_Intensity_Percentile_Value=0.9;

% Read the stored structure data
disp('Loading Data in Struct');
filepath=pwd;

% Read the initial file
filename='each_MHW_DBSCAN_Eps5_MinPts150.mat';
data_in=struct2array(load([filepath '/In/' filename]));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract parameters (extract the important parameters of the heat wave event cluster [Time, Square, Mean_Intensity] separately)
data_in_Time=data_in.Time;
data_in_Square=data_in.Square;
data_in_Mean_Intensity=data_in.Mean_Intensity;

% Calculate the total number of heat wave event clusters
MHW_step=size(data_in_Time,1);

% Sorting (sort the important parameters of the heat wave event clusters [Time, Square, Mean_Intensity])
% Sorting Time
disp('Sorting Time/Square/Mean Intensity');

% Number the heat wave event clusters in data_in_Time, [day,num] The first column is the duration of each MHW, and the second column is the MHW number
data_in_Time(:,2)=1:MHW_step;

% Sort Time
Time_sortrows=sortrows(data_in_Time, 1);

% Sorting Square
% First, generate an array Square_MAX to store the maximum area of each heat wave event cluster
Square_MAX(:,2)=1:MHW_step;
Square_MAX(:,1)=NaN;

% Calculate the maximum area of each heat wave event cluster and store it in the corresponding position in Square_MAX
for i=1:MHW_step
    nametocode_i1=['Square_MAX(' num2str(i) ',1)=max(data_in_Square.MHW_DBSCAN_Square_' num2str(i) ');'];
    eval(nametocode_i1);
end

% Sort Square
Square_sortrows=sortrows(Square_MAX, 1);

% Sorting Mean Intensity
% First, generate an array Mean_Intensity_MAX to store the maximum area of each heat wave event cluster
Mean_Intensity_MAX(:,2)=1:MHW_step;
Mean_Intensity_MAX(:,1)=NaN;

% Calculate the maximum area of each heat wave event cluster and store it in the corresponding position in Mean_Intensity_MAX
for j=1:MHW_step
    nametocode_j1=['Mean_Intensity_MAX(' num2str(j) ',1)=max(data_in_Mean_Intensity.MHW_DBSCAN_Mean_Intensity_' num2str(j) ');'];
    eval(nametocode_j1);
end

% Sort Mean_Intensity
Mean_Intensity_sortrows=sortrows(Mean_Intensity_MAX, 1);

% Store all heat wave event clusters [Time, Square_Max, Mean_Intensity_Max, Num] into Event_All
Event_All=[data_in_Time(:,1),Square_MAX(:,1),Mean_Intensity_MAX(:,1)];
Event_All(:,4)=1:MHW_step;%[Time, Square, Mean_Intensity, Num]

clear data_in_Mean_Intensity data_in_Square i j nametocode_i1 nametocode_j1
clear data_in_Time Square_MAX Mean_Intensity_MAX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KDE
disp('KDE Time/Square/Mean_Intensity');

% Set the number of kernel density estimation points
pts=MHW_step;

% Kernel density estimation for Time
[KDE_Time_Density,KDE_Time_Xi_Mesh,~]=ksdensity(Time_sortrows(:,1),'NumPoints',pts);
KDE_Time_Integral = cumtrapz(KDE_Time_Xi_Mesh, KDE_Time_Density);
KDE_Time_Result = find(KDE_Time_Integral >= KDE_Time_Percentile_Value, 1, 'first');

% Kernel density estimation for Square
[KDE_Square_Density,KDE_Square_Xi_Mesh,~]=ksdensity(Square_sortrows(:,1),'NumPoints',pts);
KDE_Square_Integral = cumtrapz(KDE_Square_Xi_Mesh, KDE_Square_Density);
KDE_Square_Result = find(KDE_Square_Integral >= KDE_Square_Percentile_Value, 1, 'first');

% Kernel density estimation for Mean_Intensity
[KDE_Mean_Intensity_Density,KDE_Mean_Intensity_Xi_Mesh,~]=ksdensity(Mean_Intensity_sortrows(:,1),'NumPoints',pts);
KDE_Mean_Intensity_Integral = cumtrapz(KDE_Mean_Intensity_Xi_Mesh, KDE_Mean_Intensity_Density);
KDE_Mean_Intensity_Result = find(KDE_Mean_Intensity_Integral >= KDE_Mean_Intensity_Percentile_Value, 1, 'first');

clear KDE_Time_Density KDE_Time_Xi_Mesh KDE_Time_Percentile_Value KDE_Time_Integral
clear KDE_Square_Density KDE_Square_Xi_Mesh KDE_Square_Percentile_Value KDE_Square_Integral
clear KDE_Mean_Intensity_Density KDE_Mean_Intensity_Xi_Mesh KDE_Mean_Intensity_Percentile_Value KDE_Mean_Intensity_Integral
clear Percentile_Value pts 

% Classification of extreme heat wave event clusters and non-extreme heat wave event clusters, classify the [Time, Square, Mean_Intensity] events into extreme and non-extreme events based on KDE
disp('Classification of Events as Normal or Extreme');

% Classification of Time
% Record the value of Time at the percentile point of the Time threshold
KDE_Time_Result_Num=find(Time_sortrows(:,2)==(KDE_Time_Result-1));
Expoint_Time_Temp=Time_sortrows(KDE_Time_Result_Num);

% Calculate the position of the next Time value greater than the threshold percentile point
for i=KDE_Time_Result_Num:MHW_step
    Expoint_Time_Temp_Next=Time_sortrows(i+1,1);
    if Expoint_Time_Temp_Next > Expoint_Time_Temp
        Expoint_Time_Lable=i+1;
        break;
    else
        Expoint_Time_Temp=Expoint_Time_Temp_Next;
    end
end

% Store events in Normal and Extreme parts in the structure Event
Event.Time.Normal=Time_sortrows(1:Expoint_Time_Lable-1,:);
Event.Time.Extreme=Time_sortrows(Expoint_Time_Lable:end,:);

% Classification of Square
% Record the value of Square at the percentile point of the Square threshold
KDE_Square_Result_Num=find(Square_sortrows(:,2)==(KDE_Square_Result));
Expoint_Square_Temp=Square_sortrows(KDE_Square_Result_Num);

% Calculate the position of the next Square value greater than the threshold percentile point
for i=KDE_Square_Result_Num:MHW_step
    Expoint_Square_Temp_Next=Square_sortrows(i+1,1);
    if Expoint_Square_Temp_Next > Expoint_Square_Temp
        Expoint_Square_Lable=i+1;
        break;
    else
        Expoint_Square_Temp=Expoint_Square_Temp_Next;
    end
end

% Store events in Normal and Extreme parts in the structure Event
Event.Square.Normal=Square_sortrows(1:Expoint_Square_Lable-1,:);
Event.Square.Extreme=Square_sortrows(Expoint_Square_Lable:end,:);

% Classification of Mean_Intensity
% Record the value of Mean_Intensity at the percentile point of the Mean_Intensity threshold
KDE_Mean_Intensity_Result_Num=find(Mean_Intensity_sortrows(:,2)==(KDE_Mean_Intensity_Result-1));
Expoint_Mean_Intensity_Temp=Mean_Intensity_sortrows(KDE_Mean_Intensity_Result_Num);

% Calculate the position of the next Mean_Intensity value greater than the threshold percentile point
for i=KDE_Mean_Intensity_Result_Num:MHW_step
    Expoint_Mean_Intensity_Temp_Next=Mean_Intensity_sortrows(i+1,1);
    if Expoint_Mean_Intensity_Temp_Next > Expoint_Mean_Intensity_Temp
        Expoint_Mean_Intensity_Lable=i+1;
        break;
    else
        Expoint_Mean_Intensity_Temp=Expoint_Mean_Intensity_Temp_Next;
    end
end

% Store events in Normal and Extreme parts in the structure Event
Event.Mean_Intensity.Normal=Mean_Intensity_sortrows(1:Expoint_Mean_Intensity_Lable-1,:);
Event.Mean_Intensity.Extreme=Mean_Intensity_sortrows(Expoint_Mean_Intensity_Lable:end,:);

clear i
clear Expoint_Time_Temp_Next Expoint_Time_Temp KDE_Time_Result Time_sortrows
clear Expoint_Square_Temp_Next Expoint_Square_Temp KDE_Square_Result Square_sortrows
clear Expoint_Mean_Intensity_Temp_Next Expoint_Mean_Intensity_Temp KDE_Mean_Intensity_Result Mean_Intensity_sortrows
clear Expoint_Time_Lable Expoint_Square_Lable Expoint_Mean_Intensity_Lable
clear KDE_Time_Result_Num KDE_Square_Result_Num KDE_Mean_Intensity_Result_Num

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding crossover events
disp('Eight Types of Crossover Events');

% 0: Normal
% 1: Extreme
% 000= Normal_Time Normal_Square Normal_Mean_Intensity
% 100= Extreme_Time Normal_Square Normal_Mean_Intensity
% 010= Normal_Time Extreme_Square Normal_Mean_Intensity
% 001= Normal_Time Normal_Square Extreme_Mean_Intensity
% 110= Extreme_Time Extreme_Square Normal_Mean_Intensity
% 101= Extreme_Time Normal_Square Extreme_Mean_Intensity
% 011= Normal_Time Extreme_Square Extreme_Mean_Intensity
% 111= Extreme_Time Extreme_Square Extreme_Mean_Intensity

% Classify into 8 categories based on the above classification logic
Event.Crossover_Event.Type_000(:,4) = setdiff(setdiff(setdiff(Event_All(:,4), Event.Time.Extreme(:,2)), Event.Square.Extreme(:,2)), Event.Mean_Intensity.Extreme(:,2));
Event.Crossover_Event.Type_100(:,4) = intersect(setdiff(Event.Time.Extreme(:,2), Event.Square.Extreme(:,2)), setdiff(Event.Time.Extreme(:,2), Event.Mean_Intensity.Extreme(:,2)));
Event.Crossover_Event.Type_010(:,4) = intersect(setdiff(Event.Square.Extreme(:,2), Event.Time.Extreme(:,2)), setdiff(Event.Square.Extreme(:,2), Event.Mean_Intensity.Extreme(:,2)));
Event.Crossover_Event.Type_001(:,4) = intersect(setdiff(Event.Mean_Intensity.Extreme(:,2), Event.Time.Extreme(:,2)), setdiff(Event.Mean_Intensity.Extreme(:,2), Event.Square.Extreme(:,2)));
Event.Crossover_Event.Type_110(:,4) = intersect(intersect(Event.Time.Extreme(:,2), Event.Square.Extreme(:,2)), setdiff(Event.Time.Extreme(:,2), Event.Mean_Intensity.Extreme(:,2)));
Event.Crossover_Event.Type_101(:,4) = intersect(intersect(Event.Time.Extreme(:,2), Event.Mean_Intensity.Extreme(:,2)), setdiff(Event.Time.Extreme(:,2), Event.Square.Extreme(:,2)));
Event.Crossover_Event.Type_011(:,4) = intersect(intersect(Event.Square.Extreme(:,2), Event.Mean_Intensity.Extreme(:,2)), setdiff(Event.Square.Extreme(:,2), Event.Time.Extreme(:,2)));
Event.Crossover_Event.Type_111(:,4) = intersect(intersect(Event.Time.Extreme(:,2), Event.Square.Extreme(:,2)), Event.Mean_Intensity.Extreme(:,2));

% Define a list of 8 event type names
TypeList = {'000' '100' '010' '001' '110' '101' '011' '111'};
sizeTypeList=size(TypeList,2);

for ii=1:sizeTypeList
    name2codeii1=char(strcat('Temp(:,4)=Event.Crossover_Event.Type_',TypeList(ii),'(:,4);'));
    eval(name2codeii1);
    sizeTemp=size(Temp,1);
    for j=1:sizeTemp
        Temp_Lable=find(Event_All(:,4) == Temp(j,4));
        Temp(j,1:3)=Event_All(Temp_Lable,1:3);
    end
    name2codeii2=char(strcat('Event.Crossover_Event.Type_',TypeList(ii),'(:,1:3)=Temp(:,1:3);'));
    eval(name2codeii2);
    clear Temp
end

clear ii j
clear Temp_Lable TypeList sizeTypeList sizeTemp name2codeii1 name2codeii2
clear MHW_step filename

%% Saving
save([filepath '/Out/KDE_Eps5_MinPts150.mat'],'data_in','Event','Event_All','Datestart','Dateend','-v7.3');