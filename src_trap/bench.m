%% This program benchmarks the density matrix solver

clear all 
clc
global G
for N = 1:1:20
for Q = 1

qr = 0.01;
omega = 0.001;
qr = qr*sqrt(omega);
Omega = 0.1;   
delta = 0.1;
varepsilon = 1;
delta_c = 1;
kappa = 1;
G = generateG(N, Q, delta, delta_c, kappa, Omega, qr, varepsilon, omega);


G(1,:)=0; % redundant equation
for m = 0:N
    n = m;
    for p = 0:Q
        q = p;
        k = m*(N+1)*(Q+1)^2+n*(Q+1)^2+p*(Q+1)+q+1;   
        G(1,k)=1;           % rho_{mm}^{p,up,p,up} to be multiplied
        G(1,k+3*(N+1)^2*(Q+1)^2)=1; % rho_{mm}^{p,dn,p,dn} to be multiplied
    end
end
RHS=zeros(4*(N+1)^2*(Q+1)^2,1);
RHS(1,1)=1; % property of the density operator, trace is one.
%% steady state solution, requiring time derivative is zero.
rho=G\RHS;
%% reshapes density operator from column to matrix form.
RMatrix_temp(:,:)=[reshape(rho(1:(N+1)^2*(Q+1)^2),(N+1)*(Q+1),(N+1)*(Q+1)).',reshape(rho(1+(N+1)^2*(Q+1)^2:2*(N+1)^2*(Q+1)^2),(N+1)*(Q+1),(N+1)*(Q+1)).';...
    reshape(rho(1+2*(N+1)^2*(Q+1)^2:3*(N+1)^2*(Q+1)^2),(N+1)*(Q+1),(N+1)*(Q+1)).',reshape(rho(1+3*(N+1)^2*(Q+1)^2:4*(N+1)^2*(Q+1)^2),(N+1)*(Q+1),(N+1)*(Q+1)).']; 
%% observables
%%%%%%%%%%%%%%%%%%%%%%%%
photonNumberMatrix = zeros(2*(N+1)*(Q+1),2*(N+1)*(Q+1)); photonSquareMatrix = photonNumberMatrix;
for m = 0:N
    for p = 0:Q
        k = m*(Q+1)+p+1; 
        photonNumberMatrix(k,k) = m;
        photonNumberMatrix(k+(N+1)*(Q+1),k+(N+1)*(Q+1)) = m;
        photonSquareMatrix(k,k) = m^2;
        photonSquareMatrix(k+(N+1)*(Q+1),k+(N+1)*(Q+1)) = m^2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%
steadystateN(N,Q)=real(trace(RMatrix_temp*photonNumberMatrix));
steadystateFluct(N,Q)=real(trace(RMatrix_temp*photonSquareMatrix))-steadystateN(N,Q)^2;
steadystateFluct(N,Q)=steadystateFluct(N,Q)/steadystateN(N,Q);% check Poissonian distribution
PartialTransposeRho(:,:)=[reshape(rho(1:(N+1)^2*(Q+1)^2),(N+1)*(Q+1),(N+1)*(Q+1)),reshape(rho(1+(N+1)^2*(Q+1)^2:2*(N+1)^2*(Q+1)^2),(N+1)*(Q+1),(N+1)*(Q+1));...
    reshape(rho(1+2*(N+1)^2*(Q+1)^2:3*(N+1)^2*(Q+1)^2),(N+1)*(Q+1),(N+1)*(Q+1)),reshape(rho(1+3*(N+1)^2*(Q+1)^2:4*(N+1)^2*(Q+1)^2),(N+1)*(Q+1),(N+1)*(Q+1))];
PTRhoEig=real(eig(PartialTransposeRho));
negativity(N,Q)=sum(abs(PTRhoEig)-PTRhoEig)/2;
clear RMatrix_temp PartialTransposeRho
N
Q
save benchmark_test_N_Q.mat
end
end


% mesh(1:8,2:9, steadystateN)
% mesh(1:8,2:9, steadystateFluct)
% mesh(1:8,2:9, negativity)
