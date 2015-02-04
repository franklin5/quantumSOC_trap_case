% clear
% clc
% clf
% figure(1)
global qr  deltac delta Omega epsilonp 
%% parameters
npara=1;
% for kappa=0.1:0.1:2
kappa=1; % chosen as energy unit
Omega=5.6; % Raman coupling strength
% for Omega=.1:0.1:6
%figure
deltac=1;
epsilonp=1;
qr=0.22; % photon recoil momentum
delta=0; % two-photon detuning
N=10; % photon number truncation
akz=-15:0.01:15;
 for kz=akz
%kz=0;    % quasimomentum of atom

steadystate2

 end
set(gca,'fontsize',16)
plot(akz,steadystateN,'r','linewidth',2)
% plot(.1:0.1:6,steadystateFluct,'linewidth',4)
% xlabel('\delta_c')
% xlabel('\Omega')
xlabel('k_z/q_r')
ylabel('<n>')
axis ([-8 8 1.5 2.5])
hold on
plot(akz,steadystateFluct./steadystateN,'--m','linewidth',2)
figure
[hAx,hLine1,hLine2] = ...
    plotyy(akz,steadystateFluct./steadystateN,...
    [aqz', aqz',aqz3', aqz4', akz'],[photon(1,:)',photon(2,:)',photon3',photon4',steadystateN']...
    ,'plot','plot');
xlabel('k_z/q_r')
ylabel(hAx(1),'fluctuation') % left y-axis
ylabel(hAx(2),'photon number') % right y-axis

% ylabel('<(\Deltan)^2>')
% ylabel('<(\Deltan)^2>/ <n>')
% figure(3)
% plot(0.1:0.1:2,steadystateFluct,'linewidth',4)
% plot(0.1:0.1:2,steadystateFluct,'linewidth',4)
% xlabel('\kappa')
% xlabel('\epsilon_p')
% xlabel('\delta_c')
% ylabel('<(\Deltan)^2> / <n>')
% saveas(figure(3),'fluc_deltac.eps','epsc')