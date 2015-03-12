%% test fluctuation parameter region with smaller than one value
%% figure wrap ups
clear
clc
clf
global qr kr eta deltac delta Omega Omega2 epsilonp 
npara=1;
kappa=1; % chosen as energy unit
Omega=5; % Raman coupling strength
Omega2=Omega;           % Omega2 is the cavity feedback. They are the same  
% for Omega=.1:0.1:6
%figure
deltac=1;
epsilonp=1.0;
eta=epsilonp/kappa;
qr=0.22; % photon recoil momentum
kr = qr;
delta=0; % two-photon detuning
N=10; % photon number truncation
akz=-50:1:-1;
photon=100*ones(4,length(akz));
 for kz=akz
steadystate2;

npara=npara+1
 end
plot(akz,steadystateFluct./steadystateN,'b',akz,ones(1,length(akz)),'r');
xlabel('k_z/q_r')
axis([min(akz) max(akz) 0.95 1.05])