function drho = recursiveC(t,rho)
global Q qr delta
G= zeros(Q,Q);
for q = 0:Q
    nq = q+1;
    G(nq,nq) = (q+0.5+delta)/1i;
    if q+1 <= Q
        G(nq,nq+1) = -qr/sqrt(2)*sqrt(q+1);
    end
    if q-1 >= 0
        G(nq,nq-1) = qr/sqrt(2)*sqrt(q);
    end
end
drho = G*rho;