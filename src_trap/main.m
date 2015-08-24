%% This program benchmarks the density matrix solver

clear 
clc
global G
npara = 1;
N = 5;
Q = 5;
% profile on
omega = 0.5;
qr = 0.22;
qr = qr*sqrt(omega);
Omega = 5;   
delta = 0.5;
varepsilon = 1;
delta_c = 1;
kappa = 1;
G = generateG(N, Q, delta, delta_c, kappa, Omega, qr, varepsilon, omega);

%% We test that the steady state evolution indeed matches the time
%% evolution results at infinite time
% rho0 = zeros(4*(N+1)^2*(Q+1)^2,1);
% rho0(1,1)=1; % initial condition of the state
% maxT = 5;
% [TimeRho, RhoT] = ode45(@timeEvoRHO, [0 maxT], rho0);
%figure(Q)
% for i = 1:length(RhoT(1,:))
%     plot(TimeRho,abs(RhoT(:,i)))
%     hold on
% end
% plot(TimeRho,abs(RhoT(:,1)),'y--')
% hold off
% figure(2)
% sump = 0;
% for p = 0:Q
%     sump = sump +  abs(RhoT(:,p*(Q+1)+p+1));
% end
% plot(TimeRho, sump)

observables
% profile viewer
% profile off