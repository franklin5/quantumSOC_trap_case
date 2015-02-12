%% semi-classical theory of photon number computation
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
    Eqplot = computeEqQuartic(qz);
    for nn=1:4
        if Eqplot(nn)~=100
            ex=Eqplot(nn);
            photon(nn,npara)=1/(1+deltac^2)*(eta^2-eta*Omega2*...
                sqrt(abs(phi_up_sqr(ex))).*sqrt(abs(phi_dn_sqr(ex))).*...
                sin(theta(ex))+(Omega2/2)^2*phi_dn_sqr(ex).*phi_up_sqr(ex));
        end
    end

%     [p,th]=pth2(Eq,qz);
%     %     display('------------------fixed point------------------')
%         [px,thx,flag,gradflag]=pthx(qz);
% 
%         hold on
%         for nn=1:length(flag)
%             orn=find(abs(exp(1i*th)-exp(1i*thx(nn)))<0.04);
%             if flag(nn)==1
%                 scatter(qz/kr,Eq(orn(1)),50,'k.') % stable point
%             else
%                 scatter(qz/kr,Eq(orn(1)),50,gradflag(nn),'v','filled') % unstable point
%             end
%         end
    