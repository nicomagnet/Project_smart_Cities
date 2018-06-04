% 5 walls: concrete (4 layesrs) + insulation (2 layers)
% glass (1 node)
% ventilation
% Inputs: outdoor temperature, solar radiation, HVAC heat flow rate
% Free running (FR) / without control

clc, clear all
% pkg load control
% Data
% ******
Sc = 7.2*6.925; Si = Sc; Sg = 10*1.01*1.48; %surface [m2]: concrete, insulation, glass
Va = Sc*2.5; %air volume[m3]

rhoa = 1.2; ca = 1000;  %indoor air density; heat capacity
Vpa = 1*Va/3600;        %infiltration and ventilation air: volume/hour

Kp = 10000; %P-controller gain: large for precision

% c: concrete; i: insulation;  g: glass
lamc = 2;       lami = 0.04;    lamg = 0.5;   %[W/m K]
rhoccc = 2.5e6; rhoici = 2020*40;   rhogcg = 1000*1000; %[J/K m3]
wc = 0.25;       wi = 0.08;      wg = 0.036;    %[m]
epswLW = 0.9;   %long wave wall emmisivity
epswSW = 0.8;   %short wave wall emmisivity

epsgLW = 0.9;   %long wave glass emmisivity
taugSW = 0.8;   %short wave glass transmittance
alphagSW = 0.2; %short wave glass absortivity

% convection coefficents
hi = 4; ho = 10;  %[W/m2 K]

% MODEL
% *****
nth = 14; nq = 16; % # of temperature node, # of flow nodes

Tm7_15=19; Tm8_9=19; Tm9_10=19; %mean temp for radiative exchange
sigma = 5.67e-8; %[W/m2 K4]
Fwg = 1/5;

% G conductance matrix
G = zeros(nq,nq);
% G(1,1)=ho*Sc; G(2,2)=lamc/(wc/8)*Sc; for i=3:9; G(i,i)=G(2,2); end;
% G(10,10)=lami/(wi/4)*Si; for i=11:13; G(i,i)=G(10,10); end; G(14,14)=hi*Si;
% G(15,15)=epswLW/(1-epswLW)*Si*4*sigma*Tm7_15^3;
% G(16,16)=Fwg*Si*4*sigma*Tm8_9^3;
% G(17,17)=epsgLW/(1-epsgLW)*Sg*sigma*Tm9_10^3;
% G(18,18)=ho*Sg; G(19,19)=lamg/(wg/2)*Sg; G(20,20)=G(19,19); G(21,21)=hi*Sg;
% G(22,22)=Vpa*rhoa*ca;

G(1,1)=ho*Sc; G(2,2)=(lamc/(wc/8))*Sc; for i=3:9; G(i,i)=G(2,2); end;
G(10,10)=(lami/(wi/4))*Si; for i=11:13; G(i,i)=G(10,10); end; G(15,15)=hi*Si;
G(14,14)=(hi*6.8*6.6);
G(16,16)=Vpa*rhoa*ca + Sg*lamg;

% C capacity matrix
C = zeros(nth);
% C(12,12)=1/4*Sc*wc*rhoccc; for i=13:15; C(i,i)=C(12,12); end;
% C(16,16)=1/2*Si*wi*rhoici; C(17,17)=C(16,16);
% C(18,18)=Sg*wg*rhogcg; C(19,19)=Va*rhoa*ca;

C(8,8)=1/4*Sc*wc*rhoccc; for i=9:11; C(i,i)=C(8,8); end;
C(12,12)=1/2*Si*wi*rhoici; C(13,13)=C(12,12);
C(14,14)=Sg*wg*rhogcg; 

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

%State-space representation
%State-space model
nnodes = size(C,1); %n° total nodes
nC = rank(C);       %n° nodes with capacity
n0 = nnodes - nC;   %n° of nodes with zero capacity

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

%Select relevant inputs and outputs
Bs = Bs(:,[[1 14 16] nq+[1 7 14]]); %inputs: [To Ti To Phiw Phii Qh]
Cs = zeros(1,nC);Cs(nC)=1;  %output
Ds = 0;

% SIMULATION
% **********

n00ReadWeatherFiles;  %read weather files

n = size(Time,1);
th = zeros(nth,n);
Qh = zeros(n,1);
TintSP = 20*ones(n,1);

%inputs
u = [Temp Temp Temp ...
  epswSW*Sc*PhiDiff taugSW*epswSW*Sg*PhiDiff ...
  Qh];

x0 = zeros(nth,n);
sys = ss(As,Bs,Cs,Ds);
[y, t, x] = lsim(sys, u, Time);
subplot(2,1,1)
plot(Time/3600,y,'g', Time/3600, Temp,'b'), xlabel('Time [h]')
title('Simulation lsim'), ylabel('T [C]')

Qh = Kp*(TintSP - y); Qh(1)=0; % thermal load
subplot(212), plot(Time/3600,Qh,'r')

% Integrate at each time step 0:dt:dt
% needs to iterate
dt = 3600/1; % simulation step 1h/dt = 3600s / dt
th = zeros(size(As,2),n);
for k = 1:n-1
 x0 = th(:,k);
 % x = lsode (@(x, t) ssm(x, t, As, Bs, u(k,:)'), x0, 0:dt:dt);
 [y, t, x] = lsim(sys, u(k:k+1,:), 0:dt:dt,x0);
 th(:,k+1) = x(2,:)';
 th(:,k+1) = min(TintSP(k+1),th(:,k+1)); % limit th_int to set point
 Qhvac(k+1) = Kp*(TintSP(k+1) - th(4,k+1));
 
%  th(:,k) = x(2,:)';
%  th(:,k) = min(TintSP(k),th(:,k)); % limit th_int to set point
%  Qhvac(k+1) = Kp*(TintSP(k) - th(4,k));
end
subplot(211), hold on, plot(Time/3600, th(4,:),'r'), hold off

subplot(212), hold on, plot(Time/3600,Qhvac), hold off
xlabel('Time [h]'), ylabel('Q_h_v_a_c [W]')