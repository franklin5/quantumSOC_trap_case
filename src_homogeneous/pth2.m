function [p,th]=pth2(Eq,qz)
global kr eta deltac delta Omega Omega2  
options=optimset('Display','off','TolFun',1e-12);
h1=qz^2/2+kr*qz+delta;
h2=qz^2/2-kr*qz-delta;
nx=1;ex=Eq(1);
for nex=1:length(Eq)  
    % ex=0 needs special treatment
    ex=Eq(nex);
    eqa=@(a) real( (ex-h1+1i*(Omega/2)^2*(1-a.^2)/(1-1i*deltac)).*...
        (ex-h2-1i*(Omega/2)^2*a.^2/(1+1i*deltac)))-(Omega/2*eta)^2/(1+deltac^2);
    imeqa=@(a) imag((ex-h1+1i*(Omega/2)^2*(1-a.^2)/(1-1i*deltac)).*...
        (ex-h2-1i*(Omega/2)^2*a.^2/(1+1i*deltac)));
%     flag1=eqa(0);
%     da=0.01;
%     ada=da:da:1-da;
%     for a=ada
%         flag2=eqa(a);
%         if sign(flag1)~=sign(flag2)
%             ax=fsolve(eqa,a,options);
%             if abs(imeqa(ax))<1e-3
%                 bx=sqrt(1-ax^2);
%                 thx=angle((ex-h1+Omega/2*1i*Omega2/2*bx^2/(1-1i*deltac))*...
%                             ax/bx/(Omega/2*eta/(1-1i*deltac)));
%                 p(nx) =bx^2-ax^2;
%                 th(nx)=thx;
%                 nx=nx+1;
%                 flag=1;
%             else
%                 flag=-1;
%             end
%         else
%             flag=-1;
%         end
%         flag1=flag2;
%     end
%     if flag==-1
%         if min(abs(eqa(ada)))<1e-6
%             flag=1;
%             index= min(abs(eqa(ada)))==abs(eqa(ada));
%             ax=ada(index);
%             bx=sqrt(1-ax^2);
%             p(nx)=bx^2-ax^2;
%             th(nx)=angle((ex-h1+Omega/2*1i*Omega2/2*bx^2/(1-1i*deltac))*...
%                             ax/bx/(Omega/2*eta/(1-1i*deltac)));
%             nx=nx+1;
%         end
%     end
        da=0.001;ada=da:da:1-da;  
        % plot(ada, eqa(ada)) 
        % plot(ada,abs(imeqa(ada)))
       
%         ax=fsolve(eqa,mean(ada(abs(imeqa(ada))<1e-1)),options);            
%         ax=fsolve(eqa,mean(ada(abs(imeqa(ada))<1e-4)),options);
        ax=fsolve(eqa,ada(abs(imeqa(ada))==min(abs(imeqa(ada)))),options);
        bx=sqrt(1-ax^2);
        thx=angle((ex-h1+Omega/2*1i*Omega2/2*bx^2/(1-1i*deltac))*...
                    ax/bx/(Omega/2*eta/(1-1i*deltac)));
%         if thx<0
%             thx=2*pi+thx;
%         end
        p(nx) =bx^2-ax^2;
        th(nx)=thx;
        nx=nx+1;


%     scatter(ada,eqa(ada),'b')
%     hold on
%     scatter(ada,imeqa(ada),'r')
%     hold off
%     grid on
%     nx
end
