%% This program benchmarks the density matrix solver

clear
clc
global G
for N = 9:2:15
for Q = 10:2:14
qr = 0.22;
Omega = 3;
delta = 1;
varepsilon = 1;
delta_c = 0.5;
kappa = 1;
G = generateG(N, Q, delta, delta_c, kappa, Omega, qr, varepsilon);
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
steadystateN(N,Q-1)=real(trace(RMatrix_temp*photonNumberMatrix));
steadystateFluct(N,Q-1)=real(trace(RMatrix_temp*photonSquareMatrix))-steadystateN(N,Q-1)^2;
steadystateFluct(N,Q-1)=steadystateFluct(N,Q-1)/steadystateN(N,Q-1);% check Poissonian distribution
PartialTransposeRho(:,:)=[reshape(rho(1:(N+1)^2*Q^2),(N+1)*Q,(N+1)*Q),reshape(rho(1+(N+1)^2*Q^2:2*(N+1)^2*Q^2),(N+1)*Q,(N+1)*Q);...
    reshape(rho(1+2*(N+1)^2*Q^2:3*(N+1)^2*Q^2),(N+1)*Q,(N+1)*Q),reshape(rho(1+3*(N+1)^2*Q^2:4*(N+1)^2*Q^2),(N+1)*Q,(N+1)*Q)];
PTRhoEig=real(eig(PartialTransposeRho));
negativity(N,Q-1)=sum(abs(PTRhoEig)-PTRhoEig)/2;
clear RMatrix_temp PartialTransposeRho
end
end
save benchmark_test_N_Q_2.mat

mesh(1:8,2:9, steadystateN)
mesh(1:8,2:9, steadystateFluct)
mesh(1:8,2:9, negativity)
