function G = generateG(N, Q, delta, delta_c, kappa, Omega, qr, varepsilon, omega)

MUU=zeros((Q+1)^2*(N+1)^2,(Q+1)^2*(N+1)^2);MUD=MUU;MDU=MUU;MDD=MUU;S0=MUU;
S1=MUU;S2=MUU;S3=MUU;S4=MUU;S5=MUU;S6=MUU;S7=MUU;S8=MUU;
for m = 0:N
    for n = 0:N
        for p = 0:Q
            for q = 0:Q
            k = m*(N+1)*(Q+1)^2+n*(Q+1)^2+p*(Q+1)+q+1;            
            %% element index is the same for k and kt
            kt = k;
            MUU(k,kt) = ((p+0.5)*omega+delta)/1i-((q+0.5)*omega+delta)/1i+1i*delta_c*(m-n)-kappa*(m+n);
            MUD(k,kt) = ((p+0.5)*omega+delta)/1i-((q+0.5)*omega-delta)/1i+1i*delta_c*(m-n)-kappa*(m+n);
            MDU(k,kt) = ((p+0.5)*omega-delta)/1i-((q+0.5)*omega+delta)/1i+1i*delta_c*(m-n)-kappa*(m+n);
            MDD(k,kt) = ((p+0.5)*omega-delta)/1i-((q+0.5)*omega-delta)/1i+1i*delta_c*(m-n)-kappa*(m+n);
            
            %% oscillator index is off by +1 or -1
            mt = m; nt = n; pt = p+1; qt = q;
            kt = mt*(N+1)*(Q+1)^2+nt*(Q+1)^2+pt*(Q+1)+qt+1;  
            if pt <= Q  
                MUU(k,kt) = -1i*qr/sqrt(2)*sqrt(p+1)/1i;
                MUD(k,kt) = -1i*qr/sqrt(2)*sqrt(p+1)/1i;
                MDU(k,kt) = +1i*qr/sqrt(2)*sqrt(p+1)/1i;
                MDD(k,kt) = +1i*qr/sqrt(2)*sqrt(p+1)/1i;
                
            end
            
            %% oscillator index is off by +1 or -1
            mt = m; nt = n; pt = p-1; qt = q;
            kt = mt*(N+1)*(Q+1)^2+nt*(Q+1)^2+pt*(Q+1)+qt+1;  
            if pt >= 0  
                MUU(k,kt) = +1i*qr/sqrt(2)*sqrt(p)/1i;
                MUD(k,kt) = +1i*qr/sqrt(2)*sqrt(p)/1i;
                MDU(k,kt) = -1i*qr/sqrt(2)*sqrt(p)/1i;
                MDD(k,kt) = -1i*qr/sqrt(2)*sqrt(p)/1i;            
            end
            
            %% oscillator index is off by +1 or -1
            mt = m; nt = n; pt = p; qt = q+1;
            kt = mt*(N+1)*(Q+1)^2+nt*(Q+1)^2+pt*(Q+1)+qt+1;  
            if qt <= Q  
                MUU(k,kt) = -1i*qr/sqrt(2)*sqrt(q+1)/1i;
                MUD(k,kt) = +1i*qr/sqrt(2)*sqrt(q+1)/1i;
                MDU(k,kt) = -1i*qr/sqrt(2)*sqrt(q+1)/1i;
                MDD(k,kt) = +1i*qr/sqrt(2)*sqrt(q+1)/1i;
            end
            
            %% oscillator index is off by +1 or -1
            mt = m; nt = n; pt = p; qt = q-1;
            kt = mt*(N+1)*(Q+1)^2+nt*(Q+1)^2+pt*(Q+1)+qt+1;  
            if qt >= 0  
                MUU(k,kt) = +1i*qr/sqrt(2)*sqrt(q)/1i;
                MUD(k,kt) = -1i*qr/sqrt(2)*sqrt(q)/1i;
                MDU(k,kt) = +1i*qr/sqrt(2)*sqrt(q)/1i;
                MDD(k,kt) = -1i*qr/sqrt(2)*sqrt(q)/1i;
            end
            
            %% photon number is off by +1 
            mt = m+1; nt = n+1; pt = p; qt = q;
            kt = mt*(N+1)*(Q+1)^2+nt*(Q+1)^2+pt*(Q+1)+qt+1;  
            if mt <= N  && nt <= N
                MUU(k,kt)=kappa*2*sqrt(m+1)*sqrt(n+1);
                MUD(k,kt)=kappa*2*sqrt(m+1)*sqrt(n+1);
                MDU(k,kt)=kappa*2*sqrt(m+1)*sqrt(n+1);
                MDD(k,kt)=kappa*2*sqrt(m+1)*sqrt(n+1);
            end
            
            %% photon number is off by +1 or -1
            mt = m+1; nt = n; pt = p; qt = q;
            kt = mt*(N+1)*(Q+1)^2+nt*(Q+1)^2+pt*(Q+1)+qt+1;  
            if mt <= N  
                MUU(k,kt) = -varepsilon*sqrt(m+1);
                MUD(k,kt) = -varepsilon*sqrt(m+1);
                MDU(k,kt) = -varepsilon*sqrt(m+1);
                MDD(k,kt) = -varepsilon*sqrt(m+1);
                S2(k,kt)  = Omega/2*sqrt(m+1)/1i;
                S4(k,kt)  = Omega/2*sqrt(m+1)/1i;
            end
            
            %% photon number is off by +1 and -1
            mt = m-1; nt = n; pt = p; qt = q;
            kt = mt*(N+1)*(Q+1)^2+nt*(Q+1)^2+pt*(Q+1)+qt+1;  
            if mt >= 0  
                MUU(k,kt) = varepsilon*sqrt(m);
                MUD(k,kt) = varepsilon*sqrt(m);
                MDU(k,kt) = varepsilon*sqrt(m);
                MDD(k,kt) = varepsilon*sqrt(m);
                S5(k,kt)  = Omega/2*sqrt(m)/1i;
                S7(k,kt)  = Omega/2*sqrt(m)/1i;
            end
            
            %% photon number is off by +1 and -1
            mt = m; nt = n+1; pt = p; qt = q;
            kt = mt*(N+1)*(Q+1)^2+nt*(Q+1)^2+pt*(Q+1)+qt+1;  
            if nt <= N  
                MUU(k,kt) = -varepsilon*sqrt(n+1);
                MUD(k,kt) = -varepsilon*sqrt(n+1);
                MDU(k,kt) = -varepsilon*sqrt(n+1);
                MDD(k,kt) = -varepsilon*sqrt(n+1);
                S1(k,kt)  = -Omega/2*sqrt(n+1)/1i;
                S6(k,kt)  = -Omega/2*sqrt(n+1)/1i;
            end

            %% photon number is off by +1 and -1
            mt = m; nt = n-1; pt = p; qt = q;
            kt = mt*(N+1)*(Q+1)^2+nt*(Q+1)^2+pt*(Q+1)+qt+1;  
            if nt >= 0  
                MUU(k,kt) = varepsilon*sqrt(n);
                MUD(k,kt) = varepsilon*sqrt(n);
                MDU(k,kt) = varepsilon*sqrt(n);
                MDD(k,kt) = varepsilon*sqrt(n);
                S3(k,kt)  = -Omega/2*sqrt(n)/1i;
                S8(k,kt)  = -Omega/2*sqrt(n)/1i;
            end            
            
            end
        end
    end
end
G=[MUU,S1,S2,S0;S3,MUD,S0,S4;S5,S0,MDU,S6;S0,S7,S8,MDD];
end