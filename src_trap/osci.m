%% This program benchmarks the density matrix solver

clear
clc
global G
N = 8;
Q = 8;
qr = 0;
Omega = 5;
alist = 0.1:0.5:5;
%for npara = 1:length(alist)     
delta = 5;
varepsilon = 0;
delta_c = 0;
%kappa = alist(npara);
kappa = 0;
G = generateG(N, Q, delta, delta_c, kappa, Omega, qr, varepsilon);

%% We test that the steady state evolution indeed matches the time
%% evolution results at infinite time
rho0 = zeros(4*(N+1)^2*Q^2,1);
rho0(1,1)=1; % initial condition of the state
maxT = 20;
[TimeRho, RhoT] = ode45(@timeEvoRHO, [0 maxT], rho0);
% figure(1)
% for i = 1:length(RhoT(1,:))
%     plot(TimeRho,abs(RhoT(:,i)))
%     hold on
% end
% hold off
plot(TimeRho,abs(RhoT(:,1)))