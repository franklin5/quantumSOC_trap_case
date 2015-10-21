clear
clc
N=100
A = zeros(N,N);
for i=1:N
    if i==1
        A(i,i+1)=-1; A(i,i)=2; 
    elseif i==N
        A(i,i-1)=-1; A(i,i)=2;
    else
        A(i,i)=2; A(i,i+1)=-1; 
        A(i,i-1)=-3; % <- if I replace -3 with -1 then the solution is well behaved.
    end
end
det(A)
cond(A)
b=A*ones(N,1);
A\b

%{

sumR = 0;
for r=0:3:3
    for m = 0:N
        n=m;
            p=Q; 
                q=p;
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
                    sumR  = sumR+RMatrix_temp(ir*(N+1)*(Q+1)+m*(Q+1)+p+1,ic*(N+1)*(Q+1)+n*(Q+1)+q+1);
                    end
              
end
sumR
sumR = 0;
for r=0:3:3
    m=N
        n=m;
            for p=0:Q 
                q=p;
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
                    sumR  = sumR+RMatrix_temp(ir*(N+1)*(Q+1)+m*(Q+1)+p+1,ic*(N+1)*(Q+1)+n*(Q+1)+q+1);
                    end
              
end
sumR
cond(G)

}%