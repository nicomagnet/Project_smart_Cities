function xdot = sw03(x, t, Time, u)
% Simple wall with capacities in all temperature nodes
% Inputs: outdoor temperature, indoor convection heat flow rate (from HVAC))

% Physical properties
% *******************
Sw = 3*3;   %wall surface [m2]
Va = 3*3*3; %air volume[m3]
% concrete:1 ; insulation:2
lam1 = 2;   lam2 = 0.04; %[W/m K]
rho1c1 = 2.5e6; rho2c2 = 2e6; rho3c3 = 1.2e3;  %[J/K m3]
w1 = 0.2;   w2 = 0.08;  %[m]
x1 = 0.05;  x2 = 0.04;  %[m]
% convection coefficents
hi = 4; ho = 10;  %[W/m2 K]
% Thermal resistances
% concrete
Rc = w1/(lam1*Sw);   Cc = Sw*w1*rho1c1;
% insulation
Ri = w2/(lam2*Sw);   Ci = Sw*w2*rho2c2;
% convection
Rvi = 1/(hi*Sw); Rvo = 1/(ho*Sw);

% Thermal circuit
%****************
nth = 7; nq = 7;  % # of temperature node, # of flow nodes
% resistances
R = zeros(nq,nq);
R(1,1) = Rvo + Rc/8; 
R(2,2) = Rc/4; R(3,3)=R(2,2); R(4,4)=R(2,2);
R(5,5) = Rc/8 + Ri/4; 
R(6,6) = Ri/2; 
R(7,7) = Ri/4 + Rvi;
G = inv(R);
% capacitances
C = zeros(nth,nth);
C(1,1) = 1/4*Sw*w1*rho1c1; C(2,2)=C(1,1); C(3,3)=C(1,1); C(4,4)=C(1,1);
C(5,5) = 1/2*Sw*w2*rho2c2; C(6,6)=C(5,5);
C(7,7) = Va*rho3c3;
% adjancecy matrix
A = eye(nq+1,nth);
A = -diff(A,1,1)';

% State-space representation
B = inv(C)*[A'*G eye(nth,nth)]; % inputs u = [b; f] size(b)=nq, size(f)=nth;
B = B(:,[1 14]);        % select the 2 relevant inputs: 1->To and 14->Qh
A = inv(C)*(-A'*G*A); 
C = zeros(1,7);C(7)=1;  % output: th(7)

ut = [interp1(Time, u(1,:), t); interp1(Time, u(2,:), t)]; 
xdot = A*x + B*ut;
endfunction