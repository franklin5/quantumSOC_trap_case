%% synthetic SOC limit

%% figure wrap ups
clear
clc
close all
global qr kr eta deltac delta Omega Omega2 epsilonp 
npara=1;
kappa=1; % chosen as energy unit
Omega=0.001; % Raman coupling strength
Omega2=Omega;           % Omega2 is the cavity feedback. They are the same  
% for Omega=.1:0.1:6
%figure
deltac=1;
epsilonp=10.0;
eta=epsilonp/kappa;
qr=0.22; % photon recoil momentum
kr = qr;
delta=0; % two-photon detuning
N=10; % photon number truncation
%akz=-30:0.013:30;
akz=-1:0.11:1;
%akz=-5:0.009:5;
photon=100*ones(4,length(akz));
for qz=akz
    % display('---------------eigenvalue method----------------')      
        Eq=computeEqQuartic(qz);
   
        [p,th]=pth2(Eq,qz);
    %     display('------------------fixed point------------------')
        [px,thx,flag,gradflag]=pthx(qz);

        hold on
        for nn=1:length(flag)
            orn=find(abs(exp(1i*th)-exp(1i*thx(nn)))<0.04);
            if flag(nn)==1
                scatter(qz,Eq(orn(1)),50,'k.') % stable point
            else
                scatter(qz,Eq(orn(1)),50,gradflag(nn),'O','filled') % unstable point
            end
        end
        drawnow
        colorbar
        caxis([0 0.3])
        colormap jet
        set(gca,'fontsize',20)
        set(gca,'nextplot','replacechildren')
    npara=npara+1
end

