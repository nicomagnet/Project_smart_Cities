% 5 walls: concrete (4 layesrs) + insulation (2 layers)
% glass (1 node)
% ventilation
% Inputs: outdoor temperature, solar radiation, HVAC heat flow rate
% Feed back(FB) / control : 1) HVAC 2) enveloppe: ventilation + solar
% Non-linear controller: only heating

clc, clear all
% pkg load control

Kp = 1000;

% Data
% ******
Sc = 5*3*3; Si = Sc; Sg = 3*3; %surface [m2]: concrete, insulation, glass
Va = 3*3*3; %air volume[m3]

rhoa = 1.2; ca = 1000;  %indoor air density; heat capacity
Vpa = 1*Va/3600;        %infiltration and ventilation air: volume/hour

% c: concrete; i: insulation;  g: glass
%lamc = 2;       lami = 0.04;    lamg = 1.2;   %[W/m K]
lamc = 0.04;       lami = 0.2;    lamg = 1.2;   %[W/m K]
%rhoccc = 2.5e6; rhoici = 2e6;   rhogcg = 2e6; %[J/K m3]
rhoccc = 2e6; rhoici = 2.5e6;   rhogcg = 2e6; %[J/K m3]
%wc = 0.2;       wi = 0.08;      wg = 0.01;    %[m]
wc = 0.08;       wi = 0.2;      wg = 0.01;    %[m]
epswLW = 0.9;   %long wave wall emmisivity
epswSW = 0.8;   %short wave wall emmisivity

epsgLW = 0.9;   %long wave glass emmisivity
taugSW = 0.8;   %short wave glass transmitance
alphagSW = 0.2; %short wave glass absortivity

% convection coefficents
hi = 4; ho = 10;  %[W/m2 K]

% MODEL
% *****
% nth = 8; nq = 12; % # of temperature node, # of flow nodes
nth=14; nq= 16;
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
G(11,11)= Vpa*rhoa*ca;
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
% Inputs to keep
b = [1 8 11 12]; f = [1 3 7 8]; 

[As, Bs] = fdae2ss(A, G, C, b, f);
Cs = zeros(1,rank(C));Cs(rank(C))=1;  %output
Ds = 0;

% SIMULATION
% **********
% Read data
n00ReadWeatherFilesSolRad;  %read weather files

n = size(Time,1);
th = zeros(nth,n);
Qh = zeros(n,1);
TintSP = 22*ones(n,1);
% Inputs vector
u = [Temp Temp Temp TintSP ...
  epswSW*Sc*PhiDif taugSW*epswSW*Sg*PhiDif alphagSW*Sg*PhiDif ...
  Qh];
  
% Integrate using lsim (linear systems)
x0 = zeros(nth,n);
sys = ss(As,Bs,Cs,Ds);
[y, t, x] = lsim(sys, u, Time);
subplot(2,1,1)
plot(Time/3600,y,'r', Time/3600, Temp,'b'), 
xlabel('Time [h]'), ylabel('T [C]')
Qh = Kp*(TintSP - y); Qh(1)=0;
subplot(212), plot(Time/3600,Qh,'r')


% Integrate at each time step 0:dt:dt
% needs to itterate
dt = 3600/1; % simulation step 1h/dt = 3600s / dt
th = zeros(rank(C),n);
for k = 1:n-1
  % Predictive control: change inputs based on future values         
  x0 = th(:,k);
  % x = lsode (@(x, t) ssm(x, t, As, Bs, u(k,:)'), x0, 0:dt:dt);
  [y, t, x] = lsim(sys, u(k:k+1,:), 0:dt:dt,x0);

  th(:,k+1) = x(2,:)';
  
  if th(4,k+1) > TintSP(k+1) %only heating
    G(12,12) = 0; u(k+1,6) = 1/3*taugSW*epswSW*Sg*PhiDif(k+1);
    else
    G(12,12) = Kp; u(k+1,6) = taugSW*epswSW*Sg*PhiDif(k+1);
  end
  
  if (th(4,k+1) > Temp(k+1)) && (th(4,k+1) > TintSP(k+1))  %free cooling
    G(11,11) = 2*Va/3600*rhoa*ca;
    else
    G(11,11) = 0.5*Va/3600*rhoa*ca;
  end
  
  if th(4,k+1) > TintSP(k+1) % shading
    u(k+1,6) = 1/3*taugSW*epswSW*Sg*PhiDif(k+1);
    else
    u(k+1,6) = taugSW*epswSW*Sg*PhiDif(k+1);
  end
  
  %  
  [As, Bs] = fdae2ss(A, G, C, b, f);
  x0 = th(:,k);
  %x = lsode (@(x, t) ssm(x, t, As, Bs, u(k,:)'), x0, 0:dt:dt);
  [y, t, x] = lsim(sys, u(k:k+1,:), 0:dt:dt,x0);

  th(:,k+1) = x(2,:)';
  Qhvac(k+1) = G(12,12)*(TintSP(k+1) - th(4,k+1));
end

subplot(211), hold on, plot(Time/3600, th(4,:),'g'), hold off
xlabel('Time [h]'), ylabel('T [C]')
subplot(212), hold on, plot(Time/3600,Qhvac,'g'), hold off
xlabel('Time [h]'), ylabel('Q_h_v_a_c [W]')

