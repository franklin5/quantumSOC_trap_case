% numerical solution to the transcendental equation of eigenenergy
% final working version
% summerized, commented, updated on 05/08/13 by Lin Dong.
clear
clc
figure(1)
clf
% close all
%% parameters that are least likely to change throughout function evaluations are set by global environment.
global kr eta deltac delta Omega Omega2 
kr=0.22;kappa=1;
eta=1;          % eta=epsilon_p/kappa
deltac=1.0;       % delta_c is the pump-cavity detuning
delta=0.0;     % delta   is the two-photon  detuning
aOmega=5.6;
aqz=-15:.01:15;
Eq_data=zeros(4,length(aOmega),length(aqz));
%Framerecording=struct('cdata', [],'colormap', []);
nOmega=1;
for nOmega=1:length(aOmega)
     Omega =aOmega(nOmega);  % Omega is the ramam coupling strength. 
%Omega=6;    
Omega2=Omega;           % Omega2 is the cavity feedback. They are the same  
    Eqplot=zeros(4,length(aqz));
    for nqz=1:length(aqz)
        qz=aqz(nqz);
%         Eq=computeEq(qz);
        Eq=computeEqQuartic(qz);
        Eq_data(:,nOmega,nqz)=Eq;
        Eqplot(:,nqz)=Eq;
    end
    
%     if Omega<4*eta
%         nqz=(length(aqz)-1)/2+1;
%         [rr1,cc1]=find(Eqplot(:,nqz)==0);
%         Eqplot(rr1,nqz)=[100;100];
%     end
    %set(gca,'nextplot','replacechildren')
    for nn=1:4
        scatter(aqz/0.22,Eqplot(nn,:),50,'k.')
        hold on
    end
    set(gca,'fontsize',20)
%     axis([min(aqz) max(aqz) -2 4])
axis([min(aqz)/0.22 max(aqz)/0.22 -3 2])
    hold off
%     plot(aqz,aqz.^2/2+sqrt((kr*aqz+delta).^2+(Omega/2)^2*eta^2/(1+deltac^2)),'r',...
%         aqz,aqz.^2/2-sqrt((kr*aqz+delta).^2+(Omega/2)^2*eta^2/(1+deltac^2)),'b',...
%         aqz,aqz.^2/2+kr*aqz+delta,'-k',aqz,aqz.^2/2-kr*aqz-delta,'-g')
    xlabel('q/k_r')
    ylabel('\epsilon')
    title(['\Omega=',num2str(Omega),...
        '     \delta_c=',num2str(deltac)])
    drawnow
%    Framerecording(nOmega) = getframe(gcf);
     %saveas(figure(1),['Omega',num2str(Omega),'deltac',num2str(deltac),'.eps'],'epsc')
     save result.mat
end

% movie2avi(Framerecording, 'loop.avi', 'compression', 'None'...
%     ,'quality',99,'fps',5);
