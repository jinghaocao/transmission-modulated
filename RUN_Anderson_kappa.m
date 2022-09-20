% ANDERSON kappa
Na = 200;
% Number of randomlized phase shifts
N_rand = 20;
% Change the parameters to suit the purposes
epsr = 0;
epsk = 0.15;

% Input: a vector containing values of alpha 
% Output: (exact solution) plot band functions (real part of w) against alpha

half=floor(Na/2);
Omega = 0.35; 
L_unit = 1; % length of unit cell
L = 1; 
alp=linspace(-pi/L,pi/L,Na);
N = 3;
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

phii = [0 pi/2 pi];
T = 2*pi/Omega;

%%% Set the discretization or truncation parameters
N_multipole=3;
d_zeta=makezetadata;
k0 = 0.001;

cx = [0,0.3,0.6]'; 
c = [cx zeros(1,N)' zeros(1,N)'];

[Jdata,Hdata] = makeJHdata0(k0,R,N_multipole);
JHijdata = makeJHijexpdata(k0,c,N_multipole);
for i=1:length(alp)
    C = makeC_1D(k0,R,alp(i),L_unit,d_zeta,Jdata,Hdata,JHijdata,N,N_multipole);
    w(i,:)=w_alp(epsr, epsk,Omega,phii, N, delta,kappa_b,rho_b,vol,C);
end
w1=real(w);
w2=imag(w);
% colorarry = {'b','r','k'};
% colorarry2 = {'b','k','r'};

figure
hold on
% Plot band functions
for j=N+1:2*N
    plot(alp,w1(:,j),'.b')
    plot(alp,w2(:,j),'.b')
    
end

% Randomnized
alp_1 = 1.153;
alp_2 = 1.37;
% alp_local = [-alp_2 -alp_1 alp_1 alp_2];
interval_1 = [alp_1-0.3:0.05:alp_1+0.3];
interval_2 = [-alp_1-0.3:0.05:-alp_1+0.3];
w_random = zeros(length(interval_1),2*N,N_rand);
w_random2 = w_random;
hold on
for i = 1:N_rand
    phii = [0+randn(1)*pi pi/2+randn(1)*pi pi+randn(1)*pi];
    for a_ind = 1:length(interval_1)
        a = interval_1(a_ind);
        C = makeC_1D(k0,R,a,L_unit,d_zeta,Jdata,Hdata,JHijdata,N,N_multipole);
        b = interval_2(a_ind);
        C_2 = makeC_1D(k0,R,b,L_unit,d_zeta,Jdata,Hdata,JHijdata,N,N_multipole);        
        w_random(a_ind,:,i)=w_alp(epsr,epsk,Omega,phii, N, delta,kappa_b,rho_b,vol,C);
        w_random2(a_ind,:,i)=w_alp(epsr,epsk,Omega,phii, N, delta,kappa_b,rho_b,vol,C_2);
    end
    plot(interval_1,real(w_random(:,5,i)))
    plot(interval_1,real(w_random(:,6,i)))
    plot(interval_1,imag(w_random(:,5,i)))
    plot(interval_1,imag(w_random(:,6,i)))
    plot(interval_2,real(w_random2(:,5,i)))
    plot(interval_2,real(w_random2(:,6,i)))
    plot(interval_2,imag(w_random2(:,5,i)))
    plot(interval_2,imag(w_random2(:,6,i)))
end

% Choose a legend and adjust the x,y lim
title("$\Omega=$ "+Omega+", $\epsilon_r=$ "+epsr+", $\epsilon_k=$ "+epsk, 'interpreter','latex')
ylabel("Re($\omega$) and Im($\omega$)",'interpreter','latex')
xlabel("\alpha")
x = [-pi/L, 0, pi/L];
y= ["$-\frac{\pi}{L}$"; "0";"$\frac{\pi}{L}$"];
xlim([-pi/L,pi/L]);
ylim([-max(max(w2))*1.2,Omega/2]);
% ylim([0,Omega/2]);

% zoom in only for kgap
% xlim([1.09,1.23]);
% ylim([-max(max(w2))*0.4,max(max(w2))*0.4])
set(gca,'XTick',x); % Change x-axis ticks
set(gca,'XTickLabel',y'); 
set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 14.5);
%axis square
pbaspect([1.2 1 1])

function w= w_alp(epsr, epsk,Omega,phii, N, delta,kappa_b,rho_b,vol,C)
    w3 = @(t) [Omega^2/4*(1+((epsk.^2-1)./(1+epsk.*cos(Omega*t+phii)).^2))];
    rhot = @(t) 1./(1+epsr.*cos(Omega*t+phii));
    sqrtkappat = @(t) [1./sqrt(1+cos(Omega*t+phii).*epsk)];
    T = 2*pi/Omega;
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    [w, V] = hill_exp(T,M,N); % what is cnd?
    w_real=real(w);
    % w_real is a vector of size 2N=6
    % Sort w_real in order to yield the correct bands
    [w_real,ind]=sort(w_real);
    w = w(ind);
end   