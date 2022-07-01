%% Numerically compute the solutions of a NxN system of Hill equations
%  \Psi''(t) + M(t)\Psi(t) = 0. M defines the matrix-valued coefficient with period T. 

function sol = hill_sol(M,M0,N,T_span,Omega)

I = eye(N,N);
Z = zeros(N,N);


T = 2*pi/Omega;
MM = @(t,y) [Z, I; -M(t), Z]*y;
MM0 = @(t,y) [Z, I; -M0(t), Z]*y;

% Thea's initial condition
% II = eye(2*N,2*N);
% X = zeros(2*N,2*N);
% for j = 1 : 2*N
%    [~,y] = ode45(MM0,[0,T],II(:,j)); 
%    X(:,j) = y(end,:)'; 
% end
% [V,d] = eig(X,'vector');
% f_exp = log(d)/1i/T;
% [f_exp_real, ind] = sort(real(f_exp),'descend');
% f_exp = f_exp(ind);
% V = V(:, ind);
% d = d(ind);

%Jinghao's initial condition
A = [Z, I; -M0(0), Z];
[V,d] = eig(A,'vector');
[f_exp, ind] = sort(real(d/1i),'descend');
V = V(:, ind);

%normalize the eigen vector
v_1 = V(:,1)/ (sqrt(2)*V(2,1));
v_2 = V(:,2)/ (sqrt(2)*V(2,2)); 
[~,sol1] = ode45(MM,T_span,v_1);
[~,sol2] = ode45(MM,T_span,v_2);

% phse = [];
% for indext = 1: length(T_span)
%     v = (sol2(indext,:)./conj(v_2'));
%     phse = [phse;v(1)];
% end
% log(phse)./T_span';
% figure
% hold on
% plot(T_span,imag(log(phse)))
% plot(T_span,imag(log(phse))./T_span')
% plot(T_span,real(log(phse)))
% exp(-1i*f_exp(1)*T_span').*sol1

sol1 = exp(-1i*f_exp(1)*T_span').*sol1;
sol2 = exp(-1i*f_exp(2)*T_span').*sol2;

sol =[sol1 sol2];

end
