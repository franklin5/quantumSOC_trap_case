MUU=zeros((N+1)^2,(N+1)^2);MUD=MUU;MDU=MUU;MDD=MUU;S0=MUU;
S1=MUU;S2=MUU;S3=MUU;S4=MUU;S5=MUU;S6=MUU;S7=MUU;S8=MUU;
for m=0:N
    for n=0:N
        k=m*(N+1)+(n+1);
        %%
        mt=m-1;nt=n;kt=mt*(N+1)+(nt+1);
        if mt>=0        
            MUU(k,kt)=+sqrt(m)*epsilonp;
            MUD(k,kt)=+sqrt(m)*epsilonp;
            S5(k,kt)=+Omega/2/1i*sqrt(m);
            MDU(k,kt)=+sqrt(m)*epsilonp;
            S7(k,kt)=+Omega/2/1i*sqrt(m);
            MDD(k,kt)=+sqrt(m)*epsilonp;
        end
        %%
        mt=m+1;nt=n;kt=mt*(N+1)+(nt+1);
        if mt<=N        
            MUU(k,kt)=-sqrt(m+1)*epsilonp;
            S2(k,kt)=+Omega/2/1i*sqrt(m+1);
            MUD(k,kt)=-sqrt(m+1)*epsilonp;
            S4(k,kt)=+Omega/2/1i*sqrt(m+1);
            MDU(k,kt)=-sqrt(m+1)*epsilonp;
            MDD(k,kt)=-sqrt(m+1)*epsilonp;
        end
        %%
        mt=m;nt=n+1;kt=mt*(N+1)+(nt+1);
        if nt<=N        
            MUU(k,kt)=-sqrt(n+1)*epsilonp;
            S1(k,kt)=-Omega/2/1i*sqrt(n+1);
            MUD(k,kt)=-sqrt(n+1)*epsilonp;
            S6(k,kt)=-Omega/2/1i*sqrt(n+1);
            MDU(k,kt)=-sqrt(n+1)*epsilonp;
            MDD(k,kt)=-sqrt(n+1)*epsilonp;
        end 
        %%
        mt=m;nt=n-1;kt=mt*(N+1)+(nt+1);
        if nt>=0        
            MUU(k,kt)=+sqrt(n)*epsilonp;
            S3(k,kt)=-Omega/2/1i*sqrt(n);
            MUD(k,kt)=+sqrt(n)*epsilonp;
            MDU(k,kt)=+sqrt(n)*epsilonp;
            S8(k,kt)=-Omega/2/1i*sqrt(n);
            MDD(k,kt)=+sqrt(n)*epsilonp;
        end  
        %%
        mt=m;nt=n;kt=mt*(N+1)+(nt+1);
        MUU(k,kt)=-deltac/1i*(m-n)-kappa*(m+n);
        MUD(k,kt)=-deltac/1i*(m-n)-kappa*(m+n)+2/1i*(qr*kz+delta);
        MDU(k,kt)=-deltac/1i*(m-n)-kappa*(m+n)-2/1i*(qr*kz+delta);
        MDD(k,kt)=-deltac/1i*(m-n)-kappa*(m+n);
        %%
        mt=m+1;nt=n+1;kt=mt*(N+1)+(nt+1);
        if mt<=N && nt<=N 
            MUU(k,kt)=kappa*2*sqrt(m+1)*sqrt(n+1);
            MUD(k,kt)=kappa*2*sqrt(m+1)*sqrt(n+1);
            MDU(k,kt)=kappa*2*sqrt(m+1)*sqrt(n+1);
            MDD(k,kt)=kappa*2*sqrt(m+1)*sqrt(n+1);
        end
    end 
end
G=[MUU,S1,S2,S0;S3,MUD,S0,S4;S5,S0,MDU,S6;S0,S7,S8,MDD];
G(1,:)=0; % redundant equation
%% trace[rho]=1
for m=0:N
    for n=0:N
        k=m*(N+1)+(n+1);
        if m==n
            G(1,k)=1;           % rho_{mn}^{up,up} to be multiplied
            G(1,k+3*(N+1)^2)=1; % rho_{mn}^{dn,dn} to be multiplied
        end
    end
end
RHS=zeros(4*(N+1)^2,1);
RHS(1,1)=1; % property of the density operator, trace is one.
%% steady state solution, requiring time derivative is zero.
rho=G\RHS;
%% reshapes density operator from column to matrix form.
RMatrix_temp(:,:)=[reshape(rho(1:(N+1)^2),N+1,N+1).',reshape(rho(1+(N+1)^2:2*(N+1)^2),N+1,N+1).';...
    reshape(rho(1+2*(N+1)^2:3*(N+1)^2),N+1,N+1).',reshape(rho(1+3*(N+1)^2:4*(N+1)^2),N+1,N+1).'];
%% observables
%%%%%%%%%%%%%%%%%%%%%%%%
photonNumberMatrix=[diag(0:N),zeros(1+N,1+N);zeros(1+N,1+N),diag(0:N)];
photonSquareMatrix=[diag(0:N).^2,zeros(1+N,1+N);zeros(1+N,1+N),diag(0:N).^2];
%%%%%%%%%%%%%%%%%%%%%%%%
MUU=zeros((N+1),(N+1));MUD=MUU;MDU=MUU;MDD=MUU;
h1=kz^2+qr*kz+delta;
h2=kz^2-qr*kz-delta;
% for m=0:N
%     i=m+1;
%     for n=0:N
%         j=n+1;
%         if m==n
%             MUU(i,j)=h1-m*deltac;
%             MDD(i,j)=h2-m*deltac;
%         elseif m==n+1
%             MUU(i,j)=1i*epsilonp*sqrt(n+1);
%             MDD(i,j)=1i*epsilonp*sqrt(n+1);
%             MDU(i,j)=Omega/2*sqrt(n+1);
%         elseif m==n-1
%             MUU(i,j)=-1i*epsilonp*sqrt(n);
%             MDD(i,j)=-1i*epsilonp*sqrt(n);
%             MUD(i,j)=Omega/2*sqrt(n);
%         end
%     end
% end
for m=0:N
    i=m+1;
    for n=0:N
        j=n+1;
        if m==n
            MUU(i,j)=h1;%-m*deltac;
            MDD(i,j)=h2;%-m*deltac;
        elseif m==n+1
%             MUU(i,j)=1i*epsilonp*sqrt(n+1);
%             MDD(i,j)=1i*epsilonp*sqrt(n+1);
            MDU(i,j)=Omega/2*sqrt(n+1);
        elseif m==n-1
%             MUU(i,j)=-1i*epsilonp*sqrt(n);
%             MDD(i,j)=-1i*epsilonp*sqrt(n);
            MUD(i,j)=Omega/2*sqrt(n);
        end
    end
end
HeffMatrix=[MUU,MUD;MDU,MDD];
steadystateH(npara)=real(trace(RMatrix_temp*HeffMatrix));
steadystateN(npara)=real(trace(RMatrix_temp*photonNumberMatrix));
steadystateFluct(npara)=real(trace(RMatrix_temp*photonSquareMatrix))-steadystateN(npara)^2;
% steadystateFluct(npara)=steadystateFluct(npara)/steadystateN(npara);% check Poissonian distribution
uu = reshape(rho(1:(N+1)^2),N+1,N+1).';
ud = reshape(rho(1+(N+1)^2:2*(N+1)^2),N+1,N+1).';
du = reshape(rho(1+2*(N+1)^2:3*(N+1)^2),N+1,N+1).';
dd = reshape(rho(1+3*(N+1)^2:4*(N+1)^2),N+1,N+1).';
rhoTA = [uu.', ud.'; du.',dd.'];
dEigrhoTA = real(eig(rhoTA));
negativity(npara) = abs(sum(dEigrhoTA(dEigrhoTA<0)));
% The following definition is only valid for pair basis state, that
% requires state A is paired with only one state B in the bipartite
% composition.
%negativity2(npara) = ...
%    sum(sum(abs(RMatrix_temp-diag(diag(RMatrix_temp)))))/2;
npara=npara+1;