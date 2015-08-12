%% This function is to compute the decomposition of block matrix indices
clear all
clc
N = 50;
Q = 50;

% m0=randi([0 N],1,1);n0=randi([0 N],1,1);
% p0=randi([0 Q],1,1);q0=randi([0 Q],1,1);
for m0 = 0:N
    for n0 = 0:N
        for p0 = 0:Q
            for q0 = 0:Q
k = m0*(N+1)*(Q+1)^2+n0*(Q+1)^2+p0*(Q+1)+q0+1;
m=floor((k-1)/((N+1)*(Q+1)^2));
n=floor((k-1-m*((N+1)*(Q+1)^2))/(Q+1)^2);
p=floor((k-1-m*((N+1)*(Q+1)^2)-n*(Q+1)^2)/(Q+1));
q=k-(m*(N+1)*(Q+1)^2+n*(Q+1)^2+p*(Q+1)+1);

if (m-m0==0 && n-n0==0 && p-p0==0 && q-q0==0) == 0
    error('error occured!!!')
    m-m0
    n-n0
    p-p0
    q-q0
end
            end
        end
    end
end

% conclusion: this is a very tricky number decomposition algorithm...
% One has to understand the difference between mod(k,N)-1 and mod(k-1,N)
% -- the latter one gives the correct answer

% clear
% clc
% N = 3;
% i0 = randi([0 N-1],1,1);
% j0 = randi([0 N-1],1,1);
% k = i0*N+j0+1;
%  
% j = mod(k-1,N);
% i = (k-(j+1))/N;
% 
% i-i0
% j-j0