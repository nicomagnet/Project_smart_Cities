%% Analysis Thermodynamic 
% 5 walls: concrete (4 layers) + insulation (1 layers)
% glass (1 node)
% Only Ventilation
% Inputs: outdoor temperature, solar radiation, HVAC heat flow rate
% Feed back(FB) / control : 1) HVAC 2) enveloppe: ventilation + solar

clc, clear all

%% Data
% Read data from

n00ReadWeatherFilesSolRad;  %read weather files
Kp = 1000; %P-controller gain: large for precision
Sc = 7.2*6.925; Si = Sc; Sg = 10*1.01*1.48; %surface [m2]
Va = Sc*2.5; %air volume[m3]

rhoa = 1.2; ca = 1000;  %indoor air density; heat capacity
Vpa = 2*Va/3600;        %infiltration and ventilation air: volume/hour

% c: concrete; i: insulation;  g: glass
lamc = 1.25;     
lami = 0.04;    lamg = 0.5;   %[W/m K]
%lamc = 4.14; rhoccc=1.2e6;
rhoccc = 6.7e5;  
rhoici = 2020*40;   rhogcg = 1000*1000; %[J/K m3]
wc = 0.25;       wi = 0.08;      wg = 0.036;    %[m]
epswLW = 0.9;   %long wave wall emmisivity
epswSW = 0.8;   %short wave wall emmisivity

epsgLW = 0.9;   %long wave glass emmisivity
taugSW = 0.8;   %short wave glass transmitance
alphagSW = 0.2; %short wave glass absortivity

% convection coefficents
hi = 4; ho = 10;  %[W/m2 K]

%% MODEL
% *****
nth = 14; nq = 16; % # of temperature node, # of flow nodes

Tm = 20;   %mean temp for radiative exchange
sigma = 5.67e-8;  %[W/m2 K4]
Fwg = 1/5;        %view factor wall - glass

% G conductance matrix
G = zeros(nq,nq);
G(1,1)=ho*Sc; G(2,2)=(lamc/(wc/8))*Sc; for i=3:9; G(i,i)=G(2,2); end;
G(10,10)=(lami/(wi/4))*Si; for i=11:13; G(i,i)=G(10,10); end; G(15,15)=hi*Si;
G(14,14)=Kp;
G(16,16)=Vpa*rhoa*ca + Sg*lamg;

% C capacity matrix
C = zeros(nth);
C(8,8)=1/4*Sc*wc*rhoccc; for i=9:11; C(i,i)=C(8,8); end;
C(12,12)=1/2*Si*wi*rhoici; C(13,13)=C(12,12);
C(14,14)=Sg*wg*rhogcg + Va*rhoa*ca; 

% A adjancy matrix
A = zeros(nq,nth);
A(1,1) = 1; A(2,1) = -1;
A(2,8) = 1; A(3,8) = -1;
A(3,2) = 1; A(4,2) = -1;
A(4,9) = 1; A(5,9) = -1;
A(5,3) = 1; A(6,3) = -1;
A(6,10) = 1; A(7,10) = -1;
A(7,4) = 1; A(8,4) = -1;
A(8,11) = 1; A(9,11) = -1;
A(9,5) = 1; A(10,5) = -1;
A(10,12) = 1; A(11,12) = -1;
A(11,6) = 1; A(12,6) = -1;
A(12,13) = 1; A(13,13) = -1;
A(13,7) = 1; A(14,14) = 1;
A(15,14) = 1; A(15,7) = -1;
A(16,14) = 1; 

% State-space representation
% State-space model
nnodes = size(C,1); %n� total nodes
nC = rank(C);       %n� nodes with capacity
n0 = nnodes - nC;   %n� of nodes with zero capacity

K = -A'*G*A;
K11 = K(1:n0,1:n0);
K12 = K(1:n0,n0+1:end);
K21 = K(n0+1:end,1:n0);
K22 = K(n0+1:end,n0+1:end);

Kb = A'*G;
Kb1 = Kb(1:n0,:);
Kb2 = Kb(n0+1:end,:);

CC = C(n0+1:end,n0+1:end);

As = inv(CC)*(-K21*inv(K11)*K12 + K22);
Bs = inv(CC)*[-K21*inv(K11)*Kb1+Kb2 -K21*inv(K11) eye(nnodes-n0,nnodes-n0)];

% Select relevant inputs and outputs
Bs = Bs(:,[[1 14 16] nq+[1 7 14]]); %inputs: [To Ti To Phiw Phii Qh]
Cs = zeros(1,nC);Cs(nC)=1;  %output
% Ds = zeros(1,length(Bs));
Ds=0;

%% SIMULATION
n = size(Time,1);
th = zeros(nth,n);
Qh = zeros(n,1);  %auxiliary sources (electrical, persons, etc.)
TintSP = 20*ones(n,1);
% Inputs
PhiTot = PhiDir + PhiDif + PhiRef;
u = [Temp TintSP Temp ...
  epswSW*Sc*PhiTot taugSW*epswSW*Sg*PhiTot ...
  Qh];

  % Integrate using lsim (linear systems)
x0 = zeros(nth,n);
sys = ss(As,Bs,Cs,Ds);
[y, t, x] = lsim(sys, u, Time);

Qh = Kp*(TintSP - y); Qh(1)=0; % thermal load

% Integrate at each time step 0:dt:dt
% needs to iterate
dt = 3600/1; % simulation step 1h/dt = 3600s / dt
th = zeros(size(As,2),n);
for k = 1:n-1
 x0 = th(:,k);
 % x = lsode (@(x, t) ssm(x, t, As, Bs, u(k,:)'), x0, 0:dt:dt);
 [y1, t, x] = lsim(sys, u(k:k+1,:), 0:dt:dt,x0);
 th(:,k+1) = x(2,:)';
 th(:,k+1) = min(TintSP(k+1),th(:,k+1)); % limit th_int to set point
 Qg(k+1) = Kp*(TintSP(k+1) - th(7,k+1));
 
%  th(:,k) = x(2,:)';
%  th(:,k) = min(TintSP(k),th(:,k)); % limit th_int to set point
%  Qhvac(k+1) = Kp*(TintSP(k) - th(7,k));
end
figure,
subplot(2,1,1)
plot(Time/3600,y,'g', Time/3600, Temp,'b',Time/3600, th(7,:),'r'), 
xlabel('Time [h]'), ylabel('T [�C]'),grid on
title('Temperature Analysis')
legend('Inside Temperature','Outside Temperature','Node 7 Temperature')
subplot(2,1,2), plot(Time/3600,Qh,'r',Time/3600,Qg),grid on
title('Energy Analysis')
legend('Sensible Heat','Auxiliary Heating')
xlabel('Time [h]'), ylabel('Q_g [W]')


