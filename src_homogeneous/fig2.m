%% fig2
clear
clc
clf
global qr kr eta deltac delta Omega Omega2 epsilonp 
kappa=1; % chosen as energy unit
%figure
deltac=1;
epsilonp=1.0;
eta=epsilonp/kappa;
qr=0.22; % photon recoil momentum
kr = qr;
delta=0; % two-photon detuning
akz = [0.001 1 6];
colorset = ['b','r','k'];
for nkz = 1:length(akz)
    kz = akz(nkz);
aOmega = 0.001:0.1:10;
npara=1;
for Omega=aOmega
    Omega2=Omega;           % Omega2 is the cavity feedback. They are the same  
    qz = kz;

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
    
    % display('---------------eigenvalue method----------------') 
    Eqplot = computeEqQuartic(qz);
    ex = mean(Eqplot(Eqplot==min(Eqplot)));
    OmegaTilde(npara) = sqrt((Omega/2)^2*1/(1+deltac^2)*(eta^2-eta*Omega2*...
                sqrt(abs(phi_up_sqr(ex))).*sqrt(abs(phi_dn_sqr(ex))).*...
                sin(theta(ex))+(Omega2/2)^2*phi_dn_sqr(ex).*phi_up_sqr(ex)));
     npara=npara+1
end
plot(aOmega, OmegaTilde,colorset(nkz),'linewidth',2)
hold on
end
% legend('k_z/q_r=0','k_z/q_r=1','k_z/q_r=10')
% xlabel('\Omega')
% ylabel('effective \Omega')
axis([0 10 0 1])