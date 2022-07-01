% Construct A0 and find alpha for which F0 has a degenerate eigenvalue
% alpha=2.39533
%%% Computes the quasifrequencies for resonator-modulated systems with
%%% a 1D linear lattice inside 2D space
% |D_i|=vol(i)
Omega = 0.3; %Omega = 0.3; gives us band gap

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
rho_0 = 9000;
rho_b = 1;
kappa_0 = 9000;
kappa_b = 1;
%%% High contrast parameter \delta
delta = rho_b/rho_0;
% weight = delta*kappa_b/rho_b/vol;

%%% Set positions of bubbles in unit cell
cx = [0,0.3,0.6]'; % column vector
c = [cx zeros(1,N)' zeros(1,N)']; % it's a matrix of dimension (N,3) with 3 columns

% Define the time modulation

epsr = 0*ones(1,N);

modr = @(epsrr,t,phi) 1./(1+epsrr*cos(Omega*t+phi));
phii = [0 pi/2 pi];
rhot = @(t) [];
for ii = 1:N
    rhot = @(t) [modr(epsr(ii),t,phii(ii))];
end
epsk= 0;
sqrtkappat = @(t) [1./sqrt(1+cos(Omega*t+phii).*epsk)];
w3 = @(t) [Omega^2/4*(1+((epsk^2-1)./(1+epsk*cos(Omega*t+phii)).^2))];
T = 2*pi/Omega;
%%% Set the discretization or truncation parameters
N_multipole=3;
%%% Precompute values of zeta function (used in computation of lattice sums, just for acceleration)
d_zeta=makezetadata;
k0 = 0.001;
[Jdata,Hdata] = makeJHdata0(k0,R,N_multipole);
JHijdata = makeJHijexpdata(k0,c,N_multipole);
%% Fold A0 and find folding numbers and eigenvalues of F0
alpha=2.39533;
C = makeC_1D(k0,R,alpha,L_unit,d_zeta,Jdata,Hdata,JHijdata,N,N_multipole);
M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
M0=M(0);
A0=[zeros(N) eye(N);-M0 zeros(N)];
% Warning: A0 has no degenerate eig, but F0 has! eig(F0) are equal to eig(A0) mod Omega, but they are all in [-Omega/2, Omega/2)

K=delta*kappa_b/rho_b;
[S,D]=eig(A0,'vector'); % S^-1*A0*S=D

m = round(imag(D)/Omega);
f0 = imag(D-1i*m*Omega);
[f0,ind]= sort(f0,'descend');
m = m(ind);
S = S(:,ind);

D=S\A0*S;
A1 = A_1_1(Omega, M0, phii, N);
A1 = S\A1*S;

A1minus = A_1_1(Omega, M0, -phii, N);
A1minus = S\A1minus*S;

F1_12=A1(1,2); F1_21=A1minus(2,1);
F1_56=A1(5,6); F1_65=A1minus(6,5);

f1=eig([0 F1_12;F1_21 0])
f2=eig([0 F1_56;F1_65 0])

% sqrt(A_1minus(5,1)*A1(1,5))
% sqrt(A_1minus(2,4)*A1(4,2))
% sqrt(A_1minus(1,5)*A1(5,1))
% sqrt(A_1minus(4,2)*A1(2,4))


%% Assemble matrix M
function out = Mfunc(t,delta,kappab,rhob,vol,C,rhot,sqrtkappat,w3)
Rho = diag(rhot(t));
Rinv = diag(1./rhot(t));
K = diag(sqrtkappat(t));
W3 = diag(w3(t));
out = delta*kappab/rhob*inv(diag(vol))*K*Rho*C*K*Rinv + W3;
end


%%

function A11=A_1_1(Omega, M0, phi, N)
    M_1_1 = zeros(N);
    for k=1:N
        for j=1:N
            if k==j
                M_1_1(k,k)=(Omega^2/2-M0(k,k))*exp(1i*phi(k))/2;
            else
                M_1_1(k,j)=-M0(k,j)*(exp(1i*phi(k))+exp(1i*phi(j)))/4;
            end
        end
    end
    A11=[zeros(N) zeros(N);-M_1_1 zeros(N)];
end





