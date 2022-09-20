% ANDERSON rho
Na = 200;
% Number of randomlized phase shifts
N_rand = 20;
% Change the parameters to suit the purposes
epsr = 0.1;
epsk = 0;

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
colorarry = {'b','r','k'};
colorarry2 = {'b','k','r'};

figure
hold on
% Plot band functions
for j=N+1:2*N
    plot(alp,w1(:,j),colorarry{j-N})
end
% Plot band gaps
for j=2*N
    a=min(w1(1:half,j));
    b=min(w1(half+1:end,j));
    line([-pi/L,0],[a,a],'Color','g');
    line([0,pi/L],[b,b],'Color','g');
end
for j=2*N-1
    a=max(w1(1:half,j));
    b=max(w1(half+1:end,j));
    line([-pi/L,0],[a,a],'Color','g');
    line([0,pi/L],[b,b],'Color','g');
end

% Randomnized
alp_1 = 1.153;
alp_2 = 1.37;
alp_local = [-alp_2 -alp_1 alp_1 alp_2];
w_random = zeros(length(alp_local),2*N,N_rand);
V_random = zeros(length(alp_local),2*N,2*N,N_rand);
hold on
for j = 1:length(alp_local)
    for i = 1:N_rand
        phii = [0+randn(1)*pi pi/2+randn(1)*pi pi+randn(1)*pi];
        C = makeC_1D(k0,R,alp_local(j),L_unit,d_zeta,Jdata,Hdata,JHijdata,N,N_multipole);
        [w_random(j,:,i),V_random(j,:,:,i)]=w_alp(epsr, epsk,Omega,phii, N, delta,kappa_b,rho_b,vol,C);
        scatter(alp_local(j),real(w_random(j,5,i)),25)
        scatter(alp_local(j),real(w_random(j,6,i)),25)
    end
    scatter(alp_local(j),real(mean(w_random(j,5,:))),200,'red','_')
    scatter(alp_local(j),real(mean(w_random(j,6,:))),200,'red','_')
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


figure
hold on
L_x = 5;
alp_show = alp_local(1);
for i = 1:N_rand
% Compute the mode for L_x extension
        % extract first unit cell
        v_1 = squeeze(V_random(1,4,1:N,i));
        % compute the extension
        v = [];
        for k = -L_x:L_x
            v = [v;v_1*exp(1i*alp_show*k)];
        end
        subplot(N_rand,1,i)
        plot([1:length(v)]-L_x*3,real(v))
end

function [w,V]= w_alp(epsr, epsk,Omega,phii, N, delta,kappa_b,rho_b,vol,C)
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
