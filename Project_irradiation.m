clear all
clc %clear console
% Load weaher data
M = dlmread('CAN_PQ_Montreal.Mirabel.csv');

% 6 is max temp
% 7 is min temp
from = 30*24;         %start time
fs=1*5*24;
period = [1:fs]';  %simulation period
Time = 3600*period;
Day = M(from+period,3); Hour = M(from+period,4);
Temp = M(from+period,7-1);      %Dry bulb temperarure [�C]
RadNDir =  M(from+period,15-1); %Direct normal solar radiation [Wh/m2]
RadHDif =  M(from+period,16-1); %Diffuse horizontal solar radiation [Wh/m2]
WDir = M(from+period,21-1);     %Wind direcrion: N=00; E=90°; S=180°; W=270�
WSpeed = M(from+period,22-1);   %wind speed [m/s]

month = M(from+period,2); %
day = M(from+period,3); %
hour = M(from+period,4); %
minute = M(from+period,5); %

% clear M

B = 90; Z = 0; L = 45; albedo = 0.2;
[PhiDir, PhiDif, PhiRef] = fSolRadTiltSurf(month, day, hour, minute, ...
  RadNDir, RadHDif, B, Z, L, albedo);
plot(Time/(24*3600), PhiDir,'b'), hold on %direct on surface
plot(Time/(24*3600), RadNDir,'g') % direct on normal to sun
plot(Time/(24*3600), PhiDif,'r')  % diffusif on surface
plot(Time/(24*3600), RadHDif,'k') % diffusif on horizontal surface
plot(Time/(24*3600), PhiRef,'m')  % reflected on surface
legend('\Phi_d_i_r','\Phi_N_d_i_r','\Phi_D_i_f', '\Phi_N_D_i_f')
title('Radiation on the surface (Winter Case)')
xlabel('Time [Days]')
ylabel('Radiation on the surface [Wh/m^{2}]')
%% plot the analysis of data 
x=Temp;
% define our filter per year
%Filtering the signal
fpy=25;
fsy=52;
l_fpy=fpy/fs;
l_fsy=fsy/fs;
l_o_fy=(l_fpy+l_fsy)/2;
delta_lambda=l_fsy-l_fpy;
Atenuation=20*log10(0.01);
M=round((3.1/delta_lambda)-1);
b = fir1(M,2*l_o_fy,'low',hanning(M+1));
yf=filter(b,1,x);
% Compare signals
figure(1);
%plot(n/fs,x,n/fs,yf,'-.r*'); grid on;legend('Original signal',...
plot(Time/(24*3600),x); grid on;
% legend('Variation Temperature');
title('Temperature in Montreal(Winter Case)')
xlabel('Time [Days]')
ylabel('Temperature[�C]')
% %% Analize the week cycle
% 
% % 1. get the information per week ( passband)
% % fundamental 52 and  ( to get all the information as is possible )
% % cut will be at 52 ( per week )
% fswl=40;
% fswh=365;
% % find M=700
% M=7000;
% delta_lambdaw=3.1/(M+1);
% l_o_lw=fswl/fs-delta_lambdaw/2;
% l_o_hw=fswh/fs-delta_lambdaw/2;
% 
% b = fir1(M,2*[l_o_lw l_o_hw],hanning(M+1));
% yf=filter(b,1,x);
% % Compare signals
% figure(2);
% %plot(n/fs,x,n/fs,yf,'-.r*'); grid on;legend('Original signal',...
% plot(x); hold on;plot(yf(round(M/2)+1:end),'r'); grid on;
% legend('Variation Temperature','Average');
% title('Temperature in Montreal - Canada')
% xlabel('Elements in hours')
% ylabel('Temperature[�C]')



