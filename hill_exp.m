%% Numerically compute the Floquet exponents of a NxN system of Hill equations
%  \Psi''(t) + M(t)\Psi(t) = 0. M defines the matrix-valued coefficient with period T. 

function [out, cnd] = hill_exp(T,M,N)
W = zeros(2*N,2*N);
Z = zeros(N,N);
I = eye(N,N);
MM = @(t,y) [Z, I; -M(t), Z]*y;
for j = 1:N
    wI0 = zeros(2*N,1); wI0(j) = 1;
    wII0 = zeros(2*N,1); wII0(N+j) = 1;
    [~, wI] = ode45(MM, [0,T], wI0);
    [~, wII] = ode45(MM, [0,T], wII0);
    W(:,j) = wI(end,:).';
    W(:,N+j) = wII(end,:).';
end

[U, D, V] = eig(W);
out = (log(diag(D))/(1i*T));
[out_real,ind] = sort(real(out),'descend');
out = out(ind);
cnd = cond(U);
end