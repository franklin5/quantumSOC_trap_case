function [px,thx,flag,gradflag]=pthx(qz)
global kr eta deltac delta Omega Omega2 
options=optimset('Display','off','TolFun',1e-12);
h1=qz^2/2+kr*qz+delta;
h2=qz^2/2-kr*qz-delta;
nx=1;
pstar=@(th) -8*(kr*qz+delta)*(1+deltac^2)/Omega^2./...
                     ((1-deltac*tan(th))./(tan(th)+deltac)+deltac);
fixedptEq=@(th) real(Omega.*sqrt(1-pstar(th).^2)-4*eta*(sin(th)+deltac*cos(th)) );
dth=0.001;
flag1=fixedptEq(0);
ath=dth:dth:pi*2-dth;
% scatter(ath,fixedptEq(ath))
% grid on                    
        for th=ath
       %for th=0.5*pi+dth:dth:pi
            flag2=fixedptEq(th);
            if sign(flag1)~=sign(flag2)  
%                 th/pi
%                 sign(flag1)
%                 sign(flag2)
                [thx(nx),fval,exitflag]=fsolve(fixedptEq,th,options);
                if exitflag~=1
                    thx(nx)=th;
                end
                %thx/pi;
                %nth
                %thx1=(asin(Omega/4/sqrt(1+deltac^2))-atan(deltac));
                %thx2=(pi-asin(Omega/4/sqrt(1+deltac^2))-atan(deltac));
                px(nx)=pstar(thx(nx));
%                 stability analysis
                %thx/pi
%                 dp=0.01;dth=0.01;
%                 pxdp=px(nx)-dp;
%                 thxdt=thx(nx)-dth;
% pertp(nx) =Omega-4*eta*(sin(thx(nx))+deltac*cos(thx(nx)))...
%     /sqrt(1-(pxdp)^2);
% pertth(nx)=h1-h2+Omega/4*px(nx)/(1+deltac^2)*...
%     (4*eta*(cos(thxdt)-deltac*sin(thxdt))...
%     /sqrt(1-px(nx)^2)+deltac*Omega);
% if pertp(nx)>0 && pertth(nx)>0
%     flag(nx)=1; % stable fixed point
% else
%     flag(nx)=-1; % unstable fixed point
% end

% stability analysis 2

            p=px(nx);
            t=thx(nx);
            f1=-p*Omega^2/2/(1+deltac^2)+...
                p*eta*Omega*deltac*cos(t)/(1+deltac^2)/sqrt(1-p^2)+...
                p*eta*Omega*sin(t)       /(1+deltac^2)/sqrt(1-p^2);
            f2=sqrt(1-p^2)/(1+deltac^2)*eta*Omega*(-cos(t)+deltac*sin(t));
            g1=Omega^2*deltac/4/(1+deltac^2)+...
                eta*Omega*(cos(t)-deltac*sin(t))/(1+deltac^2)/(1-p^2)^(3/2);
            g2=-p*eta*Omega/sqrt(1-p^2)/(1+deltac^2)*(deltac*cos(t)+sin(t));
        
%         f1=-p*Omega^2/2/(1+deltac^2)+p*eta*Omega*(deltac*cos(t)+sin(t))...
%             /(1+deltac^2)/sqrt(1-p^2);
%         f2=sqrt(1-p^2)/(1+deltac^2)*eta*Omega*(cos(t)-deltac*sin(t));
%         g1=Omega/4/(1+deltac^2)*(deltac*Omega+...
%             4*eta*(cos(t)-deltac*sin(t))/(1-p^2)^(3/2));
%         g2=-p*eta*Omega/sqrt(1-p^2)/(1+deltac^2)*(deltac*cos(t)+sin(t));
            tempeig=real(eig([f1 f2
                              g1 g2]));
            if tempeig(1)>0 || tempeig(2)>0
                flag(nx)=-1; % unstable fixed point
                gradflag(nx)=max(tempeig);
            else
                flag(nx)=+1;  % stable fixed point
                gradflag(nx)=0;
            end

%                 Omega-4*eta*(sin(thx)+deltac*cos(thx))/sqrt(1-px^2)
%                 Omega-4*eta*(sin(thx1)+deltac*cos(thx1))
%                 Omega-4*eta*(sin(thx2)+deltac*cos(thx2))
            nx=nx+1;
            end
            flag1=flag2;
            %nth=nth+1;
        end
end        
        
%           flag1=fixedptEq(0.66*pi)
%           flag2=fixedptEq(0.67*pi)
%           sign(flag1)~=sign(flag2)
%           fsolve(fixedptEq,pi*0.665)/pi