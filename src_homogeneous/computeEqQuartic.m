function Eq=computeEqQuartic(qz)
global kr eta deltac delta Omega Omega2 
Eq=100*ones(4,1);
U=Omega/2*eta/(1-1i*deltac);
V=(Omega/2)*(Omega2/2)*(-1i)/(1-1i*deltac);
h1=qz^2/2+kr*qz+delta;
h2=qz^2/2-kr*qz-delta;
W=V+conj(V);
Y=qz^4+abs(V)^2+qz^2*W;
Z=(4*qz^2+2*W);
b=-4*(h1+h2)-Z;
c=Y+(h1+h2)*Z+4*h1*h2-4*abs(U)^2;
d=-(h1+h2)*Y-h1*h2*Z+4*abs(U)^2*qz^2;
e=h1*h2*Y-abs(U)^2*qz^4;
polyEq=[4 b c d e];
temp=roots(polyEq);
ntemp=1;
for n=1:4
    if imag(temp(n))==0
        Eq(ntemp)=temp(n);
        ntemp=ntemp+1;
    end
end