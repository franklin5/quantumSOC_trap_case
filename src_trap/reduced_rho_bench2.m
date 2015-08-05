function drho = reduced_rho_bench2(t,rho)
global Q qr delta
G= zeros((Q+1)^2,(Q+1)^2);
for p = 0:Q
    np = p+1;
    for q = 0:Q
        nq = q+1;
        k = (np-1)*(Q+1)+nq;
        G(k,k) = (p+0.5+delta-(q+0.5+delta))/1i;
        if p+1 <= Q
            npt = np+1;
            kt = (npt-1)*(Q+1)+nq;
            G(k,kt) = -qr/sqrt(2)*sqrt(p+1);
        end
        if p-1 >= 0
            npt = np-1;
            kt = (npt-1)*(Q+1)+nq;
            G(k,kt) = qr/sqrt(2)*sqrt(p);
        end
        if q+1 <= Q
            nqt = nq+1;
            kt = (np-1)*(Q+1)+nqt;
            G(k,kt) = -qr/sqrt(2)*sqrt(q+1);
        end
        if q-1 >= 0
            nqt = nq-1;
            kt = (np-1)*(Q+1)+nqt;
            G(k,kt) = qr/sqrt(2)*sqrt(q);
        end
    end
end
drho = G*rho;