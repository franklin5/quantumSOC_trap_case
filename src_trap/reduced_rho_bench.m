function drho = reduced_rho_bench(t,rho)
global Omega N delta
G= zeros(4,4);
G(1,2) = -Omega/2/1i*sqrt(N);
G(1,3) = Omega/2/1i*sqrt(N);
G(2,1) = -Omega/2/1i*sqrt(N);
G(2,2) = 2*delta/1i;
G(2,4) = Omega/2/1i*sqrt(N);
G(3,1) = Omega/2/1i*sqrt(N);
G(3,3) = -2*delta/1i;
G(3,4) = -Omega/2/1i*sqrt(N);
G(4,2) = Omega/2/1i*sqrt(N);
G(4,3) = -Omega/2/1i*sqrt(N);
drho = G*rho;