% clear all
clc %clear console
% Load weaher data
% M = csvread('FRA_Lyon.074810_IWEC.csv'); % if the header is removed

% filename = 'FRA_Lyon.074810_IWEC.epw.txt';
% delimiterIn = ',';
% hederlinesIn = 8;
% M = importdata(filename,delimiterIn, hederlinesIn);
% M = M.data;

M = dlmread('CAN_PQ_Montreal.Mirabel.csv');

from = 200*24;             %start time
period = [1:31*24]';  %simulation period

Time = 3600*period;
% Day = M(from+period,3); Hour = M(from+period,4);
Temp = M(from+period,7-1); %Dry bulb temperarure [°C]
PhiDirN =  M(from+period,15-1); %Direct horizontal solar radiation [Wh/m2]
PhiDif =  M(from+period,16-1); %Diffuse solar radiation [Wh/m2]
WDir = M(from+period,21-1);   %Wind direcrion: N=00; E=90°; S=180°; W=270°
WSpeed = M(from+period,22-1); %wind speed [m/s]
PhiTot=PhiDirN+PhiDif;
clear M