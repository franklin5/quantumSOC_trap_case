%% This program benchmarks the density matrix solver

clear
clc
global G
N = 1;
Q = 2;
qr = 0.22;
Omega = 3;
alist = 0.1:0.5:5;
for npara = 1:length(alist)     
delta = 1;
varepsilon = 1;
delta_c = 0.5;
kappa = alist(npara);

G = generateG(N, Q, delta, delta_c, kappa, Omega, qr, varepsilon);

%% We test that the steady state evolution indeed matches the time
%% evolution results at infinite time
rho0 = zeros(4*(N+1)^2*Q^2,1);
rho0(1,1)=1; % initial condition of the state
maxT = 300;
[TimeRho, RhoT] = ode45(@timeEvoRHO, [0 maxT], rho0);
figure(1)
for i = 1:length(RhoT(1,:))
    plot(TimeRho,abs(RhoT(:,i)))
    hold on
end
hold off

G(1,:)=0; % redundant equation
for im = 1:N+1
    m = im-1;
    n = m;
    for p = 1:Q
        q = p;
        k = m*(N+1)*Q^2+n*Q^2+(p-1)*Q+q;   
        G(1,k)=1;           % rho_{mm}^{p,up,p,up} to be multiplied
        G(1,k+3*(N+1)^2*Q^2)=1; % rho_{mm}^{p,dn,p,dn} to be multiplied
    end
end
RHS=zeros(4*(N+1)^2*Q^2,1);
RHS(1,1)=1; % property of the density operator, trace is one.
%% steady state solution, requiring time derivative is zero.
rho=G\RHS;
figure(1)
hold on
for i = 1:length(rho)
    scatter(maxT, abs(rho(i)),'o')
end
hold off
xlabel('time')
ylabel('|\rho| components')
set(gca,'fontsize',16)
drawnow
%% reshapes density operator from column to matrix form.
RMatrix_temp(:,:)=[reshape(rho(1:(N+1)^2*Q^2),(N+1)*Q,(N+1)*Q).',reshape(rho(1+(N+1)^2*Q^2:2*(N+1)^2*Q^2),(N+1)*Q,(N+1)*Q).';...
    reshape(rho(1+2*(N+1)^2*Q^2:3*(N+1)^2*Q^2),(N+1)*Q,(N+1)*Q).',reshape(rho(1+3*(N+1)^2*Q^2:4*(N+1)^2*Q^2),(N+1)*Q,(N+1)*Q).']; 
%% observables
%%%%%%%%%%%%%%%%%%%%%%%%
photonNumberMatrix = zeros(2*(N+1)*Q,2*(N+1)*Q); photonSquareMatrix = photonNumberMatrix;
for im = 1:N+1
    m = im-1;
    for p = 1:Q
        k = m*Q+p; 
        photonNumberMatrix(k,k) = m;
        photonNumberMatrix(k+(N+1)*Q,k+(N+1)*Q) = m;
        photonSquareMatrix(k,k) = m^2;
        photonSquareMatrix(k+(N+1)*Q,k+(N+1)*Q) = m^2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%
steadystateN(npara)=real(trace(RMatrix_temp*photonNumberMatrix));
steadystateFluct(npara)=real(trace(RMatrix_temp*photonSquareMatrix))-steadystateN(npara)^2;
steadystateFluct(npara)=steadystateFluct(npara)/steadystateN(npara);% check Poissonian distribution
PartialTransposeRho(:,:)=[reshape(rho(1:(N+1)^2*Q^2),(N+1)*Q,(N+1)*Q),reshape(rho(1+(N+1)^2*Q^2:2*(N+1)^2*Q^2),(N+1)*Q,(N+1)*Q);...
    reshape(rho(1+2*(N+1)^2*Q^2:3*(N+1)^2*Q^2),(N+1)*Q,(N+1)*Q),reshape(rho(1+3*(N+1)^2*Q^2:4*(N+1)^2*Q^2),(N+1)*Q,(N+1)*Q)];
PTRhoEig=real(eig(PartialTransposeRho));
negativity(npara)=sum(abs(PTRhoEig)-PTRhoEig)/2;
end
figure(2)
plot(alist, steadystateN, 'r', alist, steadystateFluct,'b', alist, negativity, 'm')
xlabel('\kappa')
legend('<n>', 'fluctuation', 'negativity')
set(gca,'fontsize',16)
title(['cutoff number N_{photon} = ',num2str(N),',  Q_{osc} = ', num2str(Q)])