% find f0 & f1 modulating rho or kappa or both
% contrast:=delta^-1=9000; str="rho", "kappa"or "both"
% if both, phii=[phi_r; phi_k]=2x3 matrix, otherwise phi_r=phi_k=phii

alp=2.39533;
Omega=0.3; 
cx=[0 0.3 0.6]';
eps=linspace(0,0.2,300);
contrast=9000;
[A1_r,f1_r]=analytic_sol(Omega, alp,cx, phii, contrast, 'rho');
[A1_k,f1_k]=analytic_sol(Omega, alp,cx, phii, contrast, 'kappa');
[A1_b,f1_b]=analytic_sol(Omega, alp,cx, phii, contrast, 'both');

% A1_b = A1_r+A1_k
% f1_b = sqrt(f1_r^2+f1_k^2)
%% Compute f1 using analytic formula
% Input: alpha for which F0 has double degeneracy
% str= 'rho', 'kappa', or 'both'
% phii = 2xN vector [phi_r; phi_k] or 1xN vector in which case
% phi_r=phi_k=phii
% cx = center of bubbles, see RUN_resmod_1D_size.m
% contrast=rho_0/rho_b=delta^-1
function [A1,f1]=analytic_sol(Omega, alpha, cx, phii, contrast, str)

kappa_b=1; rho_b=1;
L_unit = 1; % length of unit cell?
L = 1; % I think this L is redundant, we only need R
%%% Set reciprocal vectors
b1 = pi/L_unit;
%%% Set symmetry points in reciprocal space
pointG = 0;
pointM = b1;
N = 3;
epsr = 0*ones(1,N);
%%% Radius of bubble
R=L*0.1*ones(N,1);
vol = pi*R.^2;
%%% Material parameters for bubble
delta=contrast^-1;

c = [cx zeros(1,N)' zeros(1,N)']; % it's a matrix of dimension (N,3) with 3 columns

% Define the time modulation
epsr = 0*ones(1,N);
modr = @(epsrr,t,phi) 1./(1+epsrr*cos(Omega*t+phi));
rhot = @(t) [];
for ii = 1:N
    rhot = @(t) [modr(epsr(ii),t,phii(ii))];
end

epsk= 0;
sqrtkappat = @(t) [1./sqrt(1+cos(Omega*t+phii).*epsk)];

w3 = @(t) [Omega^2/4*(1+((epsk^2-1)./(1+epsk*cos(Omega*t+phii)).^2))];
T = 2*pi/Omega;

% Fold A0 and find folding numbers and eigenvalues of F0, F1 modulating kappa only
%alpha=0.885525;
%%% Set the discretization or truncation parameters
N_multipole=3;
d_zeta=makezetadata;
k0 = 0.001;
[Jdata,Hdata] = makeJHdata0(k0,R,N_multipole);
JHijdata = makeJHijexpdata(k0,c,N_multipole);
C = makeC_1D(k0,R,alpha,L_unit,d_zeta,Jdata,Hdata,JHijdata,N,N_multipole);
M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);

M0=M(0);
A0=[zeros(N) eye(N);-M0 zeros(N)];
[m,f0]=fold(A0,Omega);
disp("m = "+m);
K=delta*kappa_b/rho_b;
[S,D]=eig(A0); % S^-1*A0*S=D
A0=diag(D); %A0 is purely imaginary, unfolded
disp("A0 = "+A0);
disp("f0 = "+f0);
w=find_double_degeneracy(f0);
a=w(1,1); b=w(1,2); % f0(a)=f0(b) is a double degeneracy
if D(a,a)*D(b,b)<0
    disp("Case 1: modulating rho is expected to yield real f1, modulating kappa is expected to yield imaginary f1") % degeneracy comes from folding lambda_i and lambda_j, both in iR_>0
else
    disp("Case 2: modulating rho is expected to yield imaginary f1, modulating kappa is expected to yield real f1") % degeneracy comes from folding lambda_i and -lambda_j,
end

[r,~]=size(phii); % 2x3 or 1x3
if r==2
    phi_r=phii(1,:); phi_k=phii(2,:);
elseif r==1
    phi_r=phii; phi_k=phii;
end

% Compute Fourier matrices of A1
    
if str=="rho"
    A1_plus1 = A_11_rho(Omega, C, K, phi_r, vol, N);
    A1_minus1 = A_11_rho(Omega, C, K, -phi_r, vol, N);
elseif str=="kappa"
    A1_plus1 = A_11_kappa(Omega, C, K, phi_k, vol, N);
    A1_minus1 = A_11_kappa(Omega, C, K, -phi_k, vol, N);
elseif str=="both"
    A1_plus1 = A_11_both(Omega, C, K, phi_r, phi_k, vol, N);
    A1_minus1 = A_11_both(Omega, C, K, -phi_r, -phi_k, vol, N);
else
    error("invalid string")
end
    
A1_plus1 = S^-1*A1_plus1*S;
A1_minus1 = S^-1*A1_minus1*S;
A1=[A1_plus1,A1_minus1];
A1=reshape(A1,[2*N,2*N,2]); % A1(:,:,1)=A1_plus1, A1(:,:,2)=A1_minus1

for i=1:2
    a=w(i,1);
    b=w(i,2); % indices of the first degenerate eigenvalue
    if abs(m(a)-m(b))==1
        F1_ab=A1(a,b,mod(m(a)-m(b),3)); % mod(m(a)-m(b),3)= 1 or 2
        F1_ba=A1(b,a,mod(m(b)-m(a),3));
        f1(:,i)=eig([0 F1_ab;F1_ba 0]);
    else 
        f1(:,i)=[0;0];
    end
    % assume f0 is purely imaginary
    disp("For the degenerate eigenvalue f0 = "+1i*imag(f0(a))+", modulating "+str+" yields f1=")
    disp(f1(:,i))
    
end

end
%% Assemble matrix M0
function out = Mfunc(t,delta,kappab,rhob,vol,C,rhot,sqrtkappat,w3)
Rho = diag(rhot(t));
Rinv = diag(1./rhot(t));
K = diag(sqrtkappat(t));
W3 = diag(w3(t));
out = delta*kappab/rhob*inv(diag(vol))*K*Rho*C*K*Rinv + W3;
end

%% find eigenvalues f of F0 and folding numbers m of eig(A0)
function [m,f]=fold(A0, Omega)
    v=eig(A0);
    im=imag(v);
    m=round(im/Omega);
    w=im-m*Omega;
    f=1i*w;
end

%% modulating kappa only
function A11=A_11_kappa(Omega, C, K, phi, vol, N)
% M_1_1:= M_1^(1)
    for i=1:N
    for j=1:N
        if i==j
            M_1_1(i,i)=(Omega^2/2-K*C(i,i)/vol(i))*exp(1i*phi(i))/2;
        else
            M_1_1(i,j)=-K*C(i,j)/vol(i)*(exp(1i*phi(i))+exp(1i*phi(j)))/4;
        end
    end
    end
    A11=[zeros(N) zeros(N);-M_1_1 zeros(N)];
end


%% modulating rho only

function A11=A_11_rho(Omega, C, K, phi_r, vol, N)
    for i=1:N
    for j=1:N
        if i==j
            M_1_1(i,j)=0;
        else
            M_1_1(i,j)=K*C(i,j)/vol(i)*(exp(1i*phi_r(i))-exp(1i*phi_r(j)))/2;
        end
    end
    end
    A11=[zeros(N) zeros(N);-M_1_1 zeros(N)];
end


%%
function A11=A_11_both(Omega, C, K, phi_r, phi_k, vol, N)
% M_1_1:= M_1^(1)
    for i=1:N
    for j=1:N
        if i==j
            M_1_1(i,i)=(Omega^2/2-K*C(i,i)/vol(i))*exp(1i*phi_k(i))/2;
        else
            M_1_1(i,j)=K*C(i,j)/vol(i)*(((exp(1i*phi_r(i))-exp(1i*phi_r(j)))/2)-(exp(1i*phi_k(i))+exp(1i*phi_k(j)))/4);
        end
    end
    end
    A11=[zeros(N) zeros(N);-M_1_1 zeros(N)];
end
% Input: (purely imaginary) vector f0
% return 2x2 matrix w containing the indices of degeneracies
% w(1,:)=indices of first degenerate eigenvalue
function w=find_double_degeneracy(f0)
w0=imag(f0);

% check if v0 has zero entries, return error if it has
for i=1:length(w0)
if w0(i)==0
    error("Entries of f0 should be purely imaginary and non zero")
end
end

N=length(f0)/2;
a=zeros(2,1);

[v,ind]=sort(w0,'descend');
for i=1:N
    if abs(v(i)-v(i+1))<=10^-4
        a(1)=i;
        a(2)=N+N-i;
        break;
    end
end
% v(a(1)) is the first degenerate eigenvalue, v(a(2)) the second
if a(1)==0
    error("Error: Check if a degeneracy exists.")
end
w=zeros(2,2);
for i=1:2
    j=1; k=1;
    while k<=2*N && j<3
        if abs(w0(k)-v(a(i)))<=10^-4
            w(i,j)=k; %record index of degenerate eig
            j=j+1; 
        end
        k=k+1;
    end
end

end