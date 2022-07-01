%%% This is used to produce Figure 3
%%% Computes w_alpha (with given alpha=2.38) in terms of epsilon

% Modulation frequency
Omega = 0.3;
% length of unit cell
L_unit = 1; 
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
rho_0 = 9000;
rho_b = 1;
kappa_0 = 9000;
kappa_b = 1;
% High contrast parameter \delta
delta = rho_b/rho_0;
% Set positions of bubbles in unit cell
cx = [0,0.3,0.6]'; 
c = [cx zeros(1,N)' zeros(1,N)']; 
% Phase shift
phii = [0 pi/2 pi];
% Period
T = 2*pi/Omega;
% Set the discretization or truncation parameters
N_multipole=3;
d_zeta=makezetadata;
k0 = 0.001;
[Jdata,Hdata] = makeJHdata0(k0,R,N_multipole);
JHijdata = makeJHijexpdata(k0,c,N_multipole);
% Define the alp_deg
alp=2.395;
% C_alpha_deg
C = makeC_1D(k0,R,alp,L_unit,d_zeta,Jdata,Hdata,JHijdata,N,N_multipole);
% Define epsilon array
eps=linspace(0,0.15,20);
% produce omega in terms of eps
w=zeros(length(eps),2*N);
for i=1:length(eps)
    % Change the inputs to suit the plot purposes
    w(i,:)=w_alp(eps(i),0,Omega,phii, N, delta,kappa_b,rho_b,vol,C);
end
% plot the 6 bands (2N=6)
figure
subplot(2,1,1)
hold on
for j=5:6
    plot(eps, real(w(:,j)))
end
legend('Re(\omega_1)','Re(\omega_2)','Location','east')
ylim([0.1185,0.1216])
xlabel('$\varepsilon_\rho$','interpreter','latex','Fontsize',15)
% % figure
% % hold on
subplot(2,1,2)
hold on
for j=5:6
    plot(eps, imag(w(:,j)))
end
%ylim for rho only
% ylim([-0.01,0.01])
%ylim for kappa only
ylim([-0.00148,0.0015])
legend('Im(\omega_1)','Im(\omega_2)','Location','northeast')
xlabel('$\varepsilon_\rho$','interpreter','latex','Fontsize',15)


function w= w_alp(epsr, epsk,Omega,phii, N, delta,kappa_b,rho_b,vol,C)
    w3 = @(t) [Omega^2/4*(1+((epsk.^2-1)./(1+epsk.*cos(Omega*t+phii)).^2))];
    rhot = @(t) 1./(1+epsr.*cos(Omega*t+phii));
    sqrtkappat = @(t) [1./sqrt(1+cos(Omega*t+phii).*epsk)];
    T = 2*pi/Omega;
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    [w, cnd] = hill_exp(T,M,N); 
    w_real=real(w);
    w_imag = imag(w);
    [w_real,ind]=sort(w_real);
    w = w(ind);
end   

%% Assemble matrix M
function out = Mfunc(t,delta,kappab,rhob,vol,C,rhot,sqrtkappat,w3)
Rho = diag(rhot(t));
Rinv = diag(1./rhot(t));
K = diag(sqrtkappat(t));
W3 = diag(w3(t));
out = delta*kappab/rhob*inv(diag(vol))*K*Rho*C*K*Rinv + W3;
end


