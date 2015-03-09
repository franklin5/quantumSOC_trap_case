% final working version of the code. 8/27/13. 
% dynamical instability.
% Lin Dong.
% in comparison with code main_homo_delta_dyn_4_2.m
% for dynamical evolutions. 
clear
clc
fign=1;
clf
global kr eta deltac delta Omega Omega2 
kr=0.22;kappa=1;
aeta=.2:.02:1;
Framerecording=struct('cdata', [],'colormap', []);
for eta=aeta          % eta=epsilon_p/kappa
    % for deltac=-1:1:1       % delta_c is the pump-cavity detuning
    deltac=1;
    delta=-0.05;
    Omega=2;   
    Omega2=Omega;           % Omega2 is the cavity feedback. They are the same  
    aqz=-4*kr:0.13/2:4*kr;

        

    for qz=aqz


    % display('---------------eigenvalue method----------------')      
        Eq=computeEqQuartic(qz);
   
        [p,th]=pth2(Eq,qz);
    %     display('------------------fixed point------------------')
        [px,thx,flag,gradflag]=pthx(qz);

        hold on
        for nn=1:length(flag)
            orn=find(abs(exp(1i*th)-exp(1i*thx(nn)))<0.04);
            if flag(nn)==1
                scatter(qz/kr,Eq(orn(1)),50,'k.') % stable point
            else
                scatter(qz/kr,Eq(orn(1)),50,gradflag(nn),'v','filled') % unstable point
            end
        end
        colorbar
        caxis([0 0.5])
        colormap jet
        set(gca,'fontsize',20)
        set(gca,'nextplot','replacechildren')
    end
    
    hold off
    axis([-4 4 -0.5 1.5])
    set(gca,'fontsize',20)
    title(['\epsilon_p=',num2str(eta,'%10.1f'),...
        ',      \Omega=',num2str(Omega),...
        ',      \delta_c=',num2str(deltac)])
    drawnow
    Framerecording(fign) = getframe(gcf);
    clf
    fign=fign+1;
end
movie2avi(Framerecording, ['stability',num2str(delta),'.avi'], 'compression', 'None'...
    ,'quality',99,'fps',5);
