%% This program benchmarks the density matrix solver
% This is to bench qr neq 0 and Omega = 0 with Chuanzhou's code.
clear
clc
global Q qr delta G
N = 1;
Q = 10;
qr = 2;
Omega = 0;     
delta = 1;
varepsilon = 0;
delta_c = 0;
kappa = 0;
G = generateG(N, Q, delta, delta_c, kappa, Omega, qr, varepsilon);

rho0 = zeros(4*(N+1)^2*(Q+1)^2,1);
rho0(1,1)=1; % initial condition of the state
maxT = 20;
[TimeRho1, RhoT] = ode45(@timeEvoRHO, [0 maxT], rho0);
reduced_rho0 = zeros((Q+1)^2,1);
reduced_rho0(1,1)=1; % initial condition of the state
[TimeRho2, reduced_RhoT] = ode45(@reduced_rho_bench2, [0 maxT], reduced_rho0);
Cq0 = zeros(Q+1,1);
Cq0(1,1) = 1; % initial condition of the state
[TimeRho3, CqT] = ode45(@recursiveC, [0 maxT], Cq0);
figure(1)
plot(TimeRho1,abs(RhoT(:,1)),'r',TimeRho2,abs(reduced_RhoT(:,1)),'b',TimeRho3, abs(CqT(:,1)).^2,'g')
legend('density matrix','reduced ODE','analytic')
