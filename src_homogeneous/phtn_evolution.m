%% photon number evolution
clear
clc
clf
global qr eta deltac delta Omega Omega2 epsilonp N kappa kz
npara=1;
kappa=1; % chosen as energy unit
Omega=3; % Raman coupling strength
Omega2=Omega;           % Omega2 is the cavity feedback. They are the same  
% for Omega=.1:0.1:6
%figure
deltac=1;
epsilonp=1.0;
eta=epsilonp/kappa;
qr=0.22; % photon recoil momentum
delta=0; % two-photon detuning
N=10; % photon number truncation
%akz=-30:0.013:30;
akz=-10:10:10;
colorset = ['b','r','k'];
%% observables
photonNumberMatrix=[diag(0:N),zeros(1+N,1+N);zeros(1+N,1+N),diag(0:N)];
rho_initial=zeros(4*(N+1)^2,1);
rho_initial(1,1)=1; % property of the density operator, trace is one. --> spin up pure state
for nkz=1:length(akz)
    kz = akz(nkz);
    [time_t,rho] = ode45(@rho_t,[0 25],rho_initial); % time evolving state from a pure spin up, and plot the time evolution of photon number
    for nt = 1:length(time_t)
        %% reshapes density operator from column to matrix form.
        rho_tmp = rho(nt,:);
        RMatrix_temp(:,:) = ...
          [reshape(rho_tmp(1:(N+1)^2),N+1,N+1).',...
           reshape(rho_tmp(1+(N+1)^2:2*(N+1)^2),N+1,N+1).';...
           reshape(rho_tmp(1+2*(N+1)^2:3*(N+1)^2),N+1,N+1).',...
           reshape(rho_tmp(1+3*(N+1)^2:4*(N+1)^2),N+1,N+1).'];
        steadystateN(nt)=real(trace(RMatrix_temp*photonNumberMatrix));
    end
    plot(time_t,steadystateN,colorset(nkz),'linewidth',2)
    hold on
    drawnow
    clear steadystateN
end
hold off


