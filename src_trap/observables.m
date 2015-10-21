G(1,:)=0; % redundant equation
for m = 0:N
    n = m;
    for p = 0:Q
        q = p;
        k = m*(N+1)*(Q+1)^2+n*(Q+1)^2+p*(Q+1)+q+1;   
        G(1,k)=1;           % rho_{mm}^{p,up,p,up} to be multiplied
        G(1,k+3*(N+1)^2*(Q+1)^2)=1; % rho_{mm}^{p,dn,p,dn} to be multiplied
    end
end
RHS=zeros(4*(N+1)^2*(Q+1)^2,1);
RHS(1,1)=1; % property of the density operator, trace is one.
% RHS = G*ones(4*(N+1)^2*(Q+1)^2,1);
%% steady state solution, requiring time derivative is zero.
rho=G\RHS
%rho = abs(rho);
% figure(1)
% hold on
% for i = 1:length(rho)
%     scatter(maxT, abs(rho(i)),'o')
% end
% hold off
% xlabel('time')
% ylabel('|\rho| components')
% set(gca,'fontsize',16)
% drawnow
%% reshapes density operator from column to matrix form.
% RMatrix_temp(:,:)=[reshape(rho(1:(N+1)^2*(Q+1)^2),(N+1)*(Q+1),(N+1)*(Q+1)).',reshape(rho(1+(N+1)^2*(Q+1)^2:2*(N+1)^2*(Q+1)^2),(N+1)*(Q+1),(N+1)*(Q+1)).';...
%    reshape(rho(1+2*(N+1)^2*(Q+1)^2:3*(N+1)^2*(Q+1)^2),(N+1)*(Q+1),(N+1)*(Q+1)).',reshape(rho(1+3*(N+1)^2*(Q+1)^2:4*(N+1)^2*(Q+1)^2),(N+1)*(Q+1),(N+1)*(Q+1)).']; 
% << __--__ The above reshape does not work as intended.
for r=0:3
    for m = 0:N
        for n = 0:N
            for p = 0:Q 
                for q = 0:Q
                    k = r*(N+1)^2*(Q+1)^2+m*(N+1)*(Q+1)^2+n*(Q+1)^2+p*(Q+1)+q+1;
                    if (r == 0) 
                        ir = 0; ic = 0;
                    elseif (r == 1)
                        ir = 0; ic = 1;
                    elseif (r == 2)
                        ir = 1; ic = 0;  
                    else
                        ir = 1; ic = 1;  
                    end
                    RMatrix_temp(ir*(N+1)*(Q+1)+m*(Q+1)+p+1,ic*(N+1)*(Q+1)+n*(Q+1)+q+1)=rho(k);
                    PartialTransposeRho(ir*(N+1)*(Q+1)+n*(Q+1)+q+1,ic*(N+1)*(Q+1)+m*(Q+1)+p+1)=rho(k);
                end
            end
        end
    end
end


%% observables
%%%%%%%%%%%%%%%%%%%%%%%%
photonNumberMatrix = zeros(2*(N+1)*(Q+1),2*(N+1)*(Q+1)); photonSquareMatrix = photonNumberMatrix;
orbitalNumberMatrix = zeros(2*(N+1)*(Q+1),2*(N+1)*(Q+1)); orbitalSquareMatrix = orbitalNumberMatrix;
nonzeros = 1;
for m = 0:N
    for p = 0:Q
        k = m*(Q+1)+p+1; 
        photonNumberMatrix(k,k) = m;
        photonNumberMatrix(k+(N+1)*(Q+1),k+(N+1)*(Q+1)) = m;
        photonSquareMatrix(k,k) = m^2;
        photonSquareMatrix(k+(N+1)*(Q+1),k+(N+1)*(Q+1)) = m^2;
        orbitalNumberMatrix(k,k) = p;
        orbitalNumberMatrix(k+(N+1)*(Q+1),k+(N+1)*(Q+1)) = p;
        orbitalSquareMatrix(k,k) = p^2;
        orbitalSquareMatrix(k+(N+1)*(Q+1),k+(N+1)*(Q+1)) = p^2;
        if (m~=0) 
            diagrhotmp(nonzeros,1)=RMatrix_temp(k,k);
            nonzeros = nonzeros+1;
            diagrhotmp(nonzeros,1)=RMatrix_temp(k+(N+1)*(Q+1),k+(N+1)*(Q+1));
            nonzeros = nonzeros+1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%
steadystateN(npara)=real(trace(RMatrix_temp*photonNumberMatrix));
steadystateFluct(npara)=real(trace(RMatrix_temp*photonSquareMatrix))-steadystateN(npara)^2;
steadystateFluct(npara)=steadystateFluct(npara)/steadystateN(npara);% check Poissonian distribution
steadystateOrbit(npara)=real(trace(RMatrix_temp*orbitalNumberMatrix));
steadystateOrbitFluc(npara)=real(trace(RMatrix_temp*orbitalSquareMatrix))-steadystateOrbit(npara)^2;
steadystateOrbitFluc(npara)=steadystateOrbitFluc(npara)/steadystateOrbit(npara);% check Poissonian distribution
% PartialTransposeRho(:,:)=[reshape(rho(1:(N+1)^2*(Q+1)^2),(N+1)*(Q+1),(N+1)*(Q+1)),reshape(rho(1+(N+1)^2*(Q+1)^2:2*(N+1)^2*(Q+1)^2),(N+1)*(Q+1),(N+1)*(Q+1));...
%     reshape(rho(1+2*(N+1)^2*(Q+1)^2:3*(N+1)^2*(Q+1)^2),(N+1)*(Q+1),(N+1)*(Q+1)),reshape(rho(1+3*(N+1)^2*(Q+1)^2:4*(N+1)^2*(Q+1)^2),(N+1)*(Q+1),(N+1)*(Q+1))];
PTRhoEig=real(eig(PartialTransposeRho));
negativity(npara)=sum(abs(PTRhoEig)-PTRhoEig)/2;
%end
% figure(2)
% plot(alist, steadystateN, 'r', alist, steadystateFluct,'b', alist, negativity, 'm')
% xlabel('\kappa')
% legend('<n>', 'fluctuation', 'negativity')
% set(gca,'fontsize',16)
steadystateN
steadystateFluct
%diag(RMatrix_temp)
steadystateOrbit
steadystateOrbitFluc
negativity
