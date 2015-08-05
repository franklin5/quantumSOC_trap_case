%% This program benchmarks the density matrix solver

clear
clc
global G Omega N delta
N = 1;
Q = 1;
qr = 0;
Omega = 2;     
delta = 8/5;
varepsilon = 0;
delta_c = 0;
kappa = 0;
G = generateG(N, Q, delta, delta_c, kappa, Omega, qr, varepsilon);

rho0 = zeros(4*(N+1)^2*Q^2,1);
rho0(1,1)=1; % initial condition of the state
maxT = 20;
[TimeRho1, RhoT] = ode45(@timeEvoRHO, [0 maxT], rho0);
reduced_rho0 = zeros(4,1);
reduced_rho0(1,1)=1; % initial condition of the state
[TimeRho2, reduced_RhoT] = ode45(@reduced_rho_bench, [0 maxT], reduced_rho0);
t = TimeRho2;
Delta = sqrt(Omega^2*N+delta^2)/2;
Cupt = cos(Delta*t)+delta/2/1i/Delta*sin(Delta*t);
Cdownt = Omega*sqrt(N)/2/1i/Delta*sin(Delta*t);
rho1 = abs(Cupt).^2;
rho2 = Cupt.*conj(Cdownt);
rho3 = conj(Cupt).*Cdownt;
rho4 = abs(Cdownt).^2;
figure(1)
plot(TimeRho1,abs(RhoT(:,1)),'r',TimeRho2,abs(reduced_RhoT(:,1)),'b',t,rho1,'g')
legend('density matrix','reduced ODE','analytic')
figure(2)
plot(TimeRho1,abs(RhoT(:,6)),'r',TimeRho2,abs(reduced_RhoT(:,2)),'b',t,abs(rho2),'g')
legend('density matrix','reduced ODE','analytic')
figure(3)
plot(TimeRho1,abs(RhoT(:,11)),'r',TimeRho2,abs(reduced_RhoT(:,3)),'b',t,abs(rho3),'g')
legend('density matrix','reduced ODE','analytic')
figure(4)
plot(TimeRho1,abs(RhoT(:,16)),'r',TimeRho2,abs(reduced_RhoT(:,4)),'b',t,rho4,'g')
legend('density matrix','reduced ODE','analytic')