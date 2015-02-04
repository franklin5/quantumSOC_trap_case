% photon number calculation
% data from result.mat
clear
clc
load result.mat
% aOmega=1:0.5:10;
for Omega =aOmega
Omega2=Omega;
nOmega=find(aOmega==Omega);
Eqplot(1:4,1:length(aqz))=Eq_data(:,nOmega,:);
photon=100*ones(4,length(aqz));
for nqz=1:length(aqz)
    qz=aqz(nqz);
    h1=qz^2/2+kr*qz+delta;
    h2=qz^2/2-kr*qz-delta;
    A=@(ex) (eta/Omega2*2)^2*...
        ((ex-h1).^2/(Omega^2/4*eta^2/(1+deltac^2))-(ex-h1)./(ex-h2));
    B=@(ex) 2*deltac*(ex-h1)/(Omega2/2)^2*Omega2/Omega;
    condition=@(ex) B(ex).^2-4*A(ex);
    sign_of_phi_dn=1; % this is questionable. 
    phi_dn_sqr=@(ex)   1/2*(B(ex)+sign_of_phi_dn*sqrt(abs(condition(ex))));
    phi_up_sqr=@(ex) (Omega/2)^2*eta^2/(1+deltac^2)*phi_dn_sqr(ex)./...
        abs(ex-h1+Omega/2*1i*Omega2/2*phi_dn_sqr(ex)/(1-1i*deltac)).^2;
    theta=@(ex) angle((ex-h1+Omega/2*1i*Omega2/2*phi_dn_sqr(ex)/(1-1i*deltac)).*...
        sqrt(abs(phi_up_sqr(ex)))./sqrt(abs(phi_dn_sqr(ex)))...
        /(Omega/2*eta/(1-1i*deltac)));
    for nn=1:4
        if Eqplot(nn,nqz)~=100
            ex=Eqplot(nn,nqz);
            photon(nn,nqz)=1/(1+deltac^2)*(eta^2-eta*Omega2*...
                sqrt(abs(phi_up_sqr(ex))).*sqrt(abs(phi_dn_sqr(ex))).*...
                sin(theta(ex))+(Omega2/2)^2*phi_dn_sqr(ex).*phi_up_sqr(ex));
        end
    end
end
figure(nOmega)
%clf
hold on
for nn=1:4
    scatter(aqz,photon(nn,:),'b.')
end
set(gca,'fontsize',16)
xlabel('q_z/k_r')
ylabel('average photon number |<c>|^2')
axis([min(aqz) max(aqz) 0 3])
title(['\Omega=',num2str(Omega),...
        '     \delta_c=',num2str(deltac)])
%saveas(figure(nOmega),['photon_',num2str(nOmega),'.eps'],'epsc')
%save photon.mat
end