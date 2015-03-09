%% figure wrap ups
clear
clc
close all
global qr kr eta deltac delta Omega Omega2 epsilonp 
npara=1;
kappa=1; % chosen as energy unit
Omega=0.1; % Raman coupling strength
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
%akz=-30:0.013:30;
akz=-10:1.1:10;
%akz=-5:0.009:5;
photon=100*ones(4,length(akz));
for kz=akz
    steadystate2;
    photon2;
    npara=npara+1
end
%   set(gca,'fontsize',16)
%  if Omega < 4*epsilonp
%      plot(akz,steadystateN,'r--', ...
%          akz, photon(1,:),'b',...
%          akz,photon(2,:),'b','linewidth',2)
%  else
%      r3 = find(photon(3,:)~=100);aqz3 = akz(r3);photon3 = photon(3,r3);
%      r4 = find(photon(4,:)~=100);aqz4 = akz(r4);photon4 = photon(4,r4);
%      plot(akz,steadystateN,'r--','linewidth',2)
%      hold on
%      scatter(akz, photon(1,:),'b.')
%      scatter(akz,photon(2,:),'b.')     
%      scatter(aqz3,photon3,'b.')
%      scatter(aqz4,photon4,'b.')
%      hold off
%  end
hold on 
 plot(akz,steadystateN,'r--','linewidth',2)
hold off
%axis([min(akz) max(akz) 0.3 0.7])
axis([min(akz) max(akz) 0 1])
xlabel('k_z/q_r')
%legend('|<c>|^2','Tr[\rho n]')
%save(['Omega_',num2str(Omega),'.mat'])
figure(2)
[hAx,hLine1,hLine2] = ...
    plotyy(akz,steadystateFluct./steadystateN,...
    akz,negativity);
xlabel('k_z/q_r')
ylabel(hAx(1),'fluctuation') % left y-axis
ylabel(hAx(2),'negativity') % right y-axis
