%%% Computes w_alpha (with given alpha=2.38) in terms of epsilon and gives
%%% polynomial approximation 
%%% yields numerical solution using multipole expansion & polyfit
%%% Recall f=iw is the Floquet exponent, f=f0+eps*f1 to first order

alp=2.39533;
Omega=0.3; cx=[0 0.3 0.6]';
eps=linspace(0,0.2,200);
phii=[0 pi/2 pi];
[f_num,f1_rho, f1_kappa, f1_both]=numerical_sol(Omega, alp,eps, cx, phii, 9000);
f1_both

%% compute f1 numerically
% cx= column vector giving centers of bubbles
function [f_num,f1_rho, f1_kappa, f1_both]=numerical_sol(Omega, alp, eps, cx, phii, contrast)
L_unit = 1; % length of unit cell
L = 1; 
%%% Set reciprocal vectors
b1 = pi/L_unit;
%%% Set symmetry points in reciprocal space
pointG = 0;
pointM = b1;
N = 3;
%%% Radius of bubble
R=L*0.1*ones(N,1);
vol = pi*R.^2;
%%% Material parameters for bubble
rho_b = 1;
kappa_b = 1;
%%% High contrast parameter \delta
delta = contrast^-1;
c = [cx zeros(1,N)' zeros(1,N)']; % it's a matrix of dimension (N,3) with 3 columns
T = 2*pi/Omega;

%%% Set the discretization or truncation parameters
N_multipole=3;
d_zeta=makezetadata;
k0 = 0.001;
[Jdata,Hdata] = makeJHdata0(k0,R,N_multipole);
JHijdata = makeJHijexpdata(k0,c,N_multipole);
C = makeC_1D(k0,R,alp,L_unit,d_zeta,Jdata,Hdata,JHijdata,N,N_multipole);

w=zeros(length(eps),2*N);

[r,~]=size(phii); % 2xN or 1xN
if r==2
    phi_r=phii(1,:); phi_k=phii(2,:);
elseif r==1
    phi_r=phii; phi_k=phii;
end

for i=1:length(eps)
    % w_alp(epsr, epsk,Omega,phii, N, delta,kappa_b,rho_b,vol,C)
    w(i,:,1)=w_alp(eps(i),0,Omega, phi_r, phi_k, N, delta,kappa_b,rho_b,vol,C); %rho
    w(i,:,2)=w_alp(0,eps(i),Omega, phi_k, phi_k, N, delta,kappa_b,rho_b,vol,C); %kappa
    w(i,:,3)=w_alp(eps(i),eps(i),Omega, phi_r, phi_k, N, delta,kappa_b,rho_b,vol,C);%both
end

for k=1:3
    w_real=[];
    w_imag=[];
    for j = 1:6
     w_real=[w_real; polyfit(eps,real(w(:,j,k)),6)];
     w_imag=[w_imag; polyfit(eps,imag(w(:,j,k)),6)];
    end
    % f=iw=iw_real-w_imag
    f_num(:,:,k)=1i*w_real-w_imag; 
end
    f0=f_num(:,7,1);
    f1_rho=f_num(:,6,1);
    f1_kappa=f_num(:,6,2);
    f1_both=f_num(:,6,3);
% polyfit(x,y,n) returns the coef of a polynomial
% p(1)*x^n+...+p(n)*x+p(n+1) that approximates y 
end
%% 
% For a given alpha, compute w_alpha given C = C_alpha, 
% sort w wrt real part in descending order 
function w = w_alp(epsr, epsk,Omega,phi_r, phi_k, N, delta,kappa_b,rho_b,vol,C)
    w3 = @(t) [Omega^2/4*(1+((epsk.^2-1)./(1+epsk.*cos(Omega*t+phi_k)).^2))];
    rhot = @(t) 1./(1+epsr.*cos(Omega*t+phi_r));
    sqrtkappat = @(t) [1./sqrt(1+cos(Omega*t+phi_k).*epsk)];
    T = 2*pi/Omega;
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    [w, cnd] = hill_exp(T,M,N); % cnd=condition number of M
    [w_real,ind]=sort(real(w),"descend");
    w=w(ind);    
end   

%% Assemble matrix M
function out = Mfunc(t,delta,kappab,rhob,vol,C,rhot,sqrtkappat,w3)
Rho = diag(rhot(t));
Rinv = diag(1./rhot(t));
K = diag(sqrtkappat(t));
W3 = diag(w3(t));
out = delta*kappab/rhob*inv(diag(vol))*K*Rho*C*K*Rinv + W3;
end


