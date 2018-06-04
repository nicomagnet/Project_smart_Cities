% 5 walls: concrete (4 layesrs) + insulation (2 layers)
% glass (1 node)
% ventilation
% Inputs: outdoor temperature, solar radiation, HVAC heat flow rate
% Feed back(FB) / control : 1) HVAC 2) enveloppe: ventilation + solar

clc, clear all
% pkg load control

Kp = 1000;

% Data
% ******
Sc = 5*3*3; Si = Sc; Sg = 3*3; %surface [m2]: concrete, insulation, glass
Va = 3*3*3; %air volume[m3]

rhoa = 1.2; ca = 1000;  %indoor air density; heat capacity
Vpa = 2*Va/3600;        %infiltration and ventilation air: volume/hour

% c: concrete; i: insulation;  g: glass
lamc = 2;       lami = 0.04;    lamg = 1.2;   %[W/m K]
rhoccc = 2.5e6; rhoici = 2e6;   rhogcg = 2e6; %[J/K m3]
wc = 0.2;       wi = 0.08;      wg = 0.01;    %[m]
epswLW = 0.9;   %long wave wall emmisivity
epswSW = 0.8;   %short wave wall emmisivity

epsgLW = 0.9;   %long wave glass emmisivity
taugSW = 0.8;   %short wave glass transmitance
alphagSW = 0.2; %short wave glass absortivity

% convection coefficents
hi = 4; ho = 10;  %[W/m2 K]

% MODEL
% *****
nth = 8; nq = 12; % # of temperature node, # of flow nodes

Tm = 20 + 273;   %mean temp for radiative exchange
sigma = 5.67e-8;  %[W/m2 K4]
Fwg = 1/6;        %view factor wall - glass

% G conductance matrix
G = zeros(nq,nq);
G(1,1)=ho*Sc; G(2,2)=lamc/(wc/2)*Sc; G(3,3) = G(2,2);
G(4,4)=lami/(wi/2)*Si; G(5,5) = G(4,4); G(6,6) = hi*Sc;
G7a = epswLW/(1-epswLW)*Si*4*sigma*Tm^3;
G7b = Fwg*Si*4*sigma*Tm^3;
G7c = epsgLW/(1-epsgLW)*Sg*sigma*Tm^3;
G(7,7) = 1/(1/G7a + 1/G7b + 1/G7c);
G(8,8) = 1/(1/(ho*Sg) + wg/(2*lamg*Sg)); G(9,9)=lamg/(wg/2)*Sg; 
G(10,10) = hi*Sg;
G(11,11)=Vpa*rhoa*ca;
G(12,12) = Kp;  %amplification of controller

% C capacity matrix
C = zeros(nth);
C(5,5) = Sc*wc*rhoccc;
C(6,6) = Si*wi*rhoici;
C(7,7) = Sg*wg*rhogcg; 
C(8,8)=Va*rhoa*ca;

% A adjancy matrix
A = zeros(nq,nth);
A(1,1)=1;
A(2,1)=-1;  A(2,5)=1;
A(3,2)=1;   A(3,5)=-1;
A(4,2)=-1;  A(4,6)=1;
A(5,3)=1;   A(5,6)=-1;
A(6,3)=-1;  A(6,8)=1;
A(7,3)=-1;  A(7,4)=1;
A(8,7)=1;
A(9,4)=1;   A(9,7)=-1;
A(10,4)=-1; A(10,8)=1;
A(11,8)=1;
A(12,8)=1;

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
Bs = Bs(:,[[1 8 11 12] nq+[1 3 7 8]]); %inputs: [To To To Phiw Phii Phig Qh]
Cs = zeros(1,nC);Cs(nC)=1;  %output

% SIMULATION
% **********
% Read data
n00ReadWeatherFiles;  %read weather files

n = size(Time,1);
th = zeros(nth,n);
Qh = zeros(n,1);
TintSP = 20*ones(n,1);
% Inputs
u = [Temp Temp Temp TintSP ...
  epswSW*Sc*PhiDif taugSW*epswSW*Sg*PhiTot alphagSW*Sg*PhiTot ...
  Qh];

 % Integrate using lsim (linear systems)
x0 = zeros(nth,n);
sys = ss(As,Bs,Cs);
[y, t, x] = lsim(sys, u, Time);
subplot(2,1,1)
plot(Time/3600,y,'g', Time/3600, Temp,'b'), 
xlabel('Time [h]'), ylabel('T [C]')
Qh = Kp*(TintSP - y); Qh(1)=0;
subplot(212), plot(Time/3600,Qh,'r')

% Integrate at each time step 0:dt:dt
% needs to itterate
dt = 3600/1; % simulation step 1h/dt = 3600s / dt
th = zeros(size(As,2),n);
for k = 1:n-1
 x0 = th(:,k);
 x = lsode (@(x, t) ssm(x, t, As, Bs, u(k,:)'), x0, 0:dt:dt);
 th(:,k+1) = x(2,:)';
 Qhvac(k+1) = Kp*(TintSP(k+1) - th(4,k+1));
end
subplot(211), hold on, plot(Time/3600, th(4,:),'r'), hold off

subplot(212), hold on, plot(Time/3600,Qhvac), hold off
xlabel('Time [h]'), ylabel('Q_h_v_a_c [W]')
