% Input: a vector containing values of alpha 
% Output: (exact solution) plot band functions (real part of w) against alpha

function w=plot_bands(Na,epsr,epsk)

Omega = 0.4; %Omega = 0.3; gives us band gap

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
% Omega=0.4;
% cx = [0,0.26,0.52]';
c = [cx zeros(1,N)' zeros(1,N)'];

[Jdata,Hdata] = makeJHdata0(k0,R,N_multipole);
JHijdata = makeJHijexpdata(k0,c,N_multipole);

for i=1:length(alp)
    C = makeC_1D(k0,R,alp(i),L_unit,d_zeta,Jdata,Hdata,JHijdata,N,N_multipole);
    w(i,:)=w_alp(epsr, epsk,Omega,phii, N, delta,kappa_b,rho_b,vol,C);
end

w1=real(w);
w2=imag(w);

figure;
hold on
for j=N+1:2*N
    plot(alp,w1(:,j),'b')
end
plot(2.3838,0.12,'*','MarkerSize',10)
lgd = legend('Band fuctions','','','\alpha_{deg}')
lgd.FontSize = 10;
legend('Location','southeast');
% title("$\Omega=$ "+Omega+", $\epsilon_r=$ "+epsr+", $\epsilon_k=$ "+epsk, 'interpreter','latex')
ylabel("Re($\omega$)",'interpreter','latex')
xlabel("\alpha")

x = [-pi/L, 0, pi/L];
y= ["$-\frac{\pi}{L}$";  "0"; "$\frac{\pi}{L}$"];

xlim([-pi/L,pi/L]);
ylim([0,Omega/2]);
set(gca,'XTick',x); % Change x-axis ticks
set(gca,'XTickLabel',y'); 
set(gca,'TickLabelInterpreter','latex')
set(gca, 'FontSize', 14.5);
%axis square
pbaspect([1.2 1 1])
end

function w_real= w_alp(epsr, epsk,Omega,phii, N, delta,kappa_b,rho_b,vol,C)
    w3 = @(t) [Omega^2/4*(1+((epsk.^2-1)./(1+epsk.*cos(Omega*t+phii)).^2))];
    rhot = @(t) 1./(1+epsr.*cos(Omega*t+phii));
    sqrtkappat = @(t) [1./sqrt(1+cos(Omega*t+phii).*epsk)];
    T = 2*pi/Omega;
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    [w_real, cnd] = hill_exp(T,M,N); % what is cnd?
    w_real=real(w_real);
    % w_real is a vector of size 2N=6
    % Sort w_real in order to yield the correct bands
    w_real=sort(w_real);
end   

%% Assemble matrix M
function out = Mfunc(t,delta,kappab,rhob,vol,C,rhot,sqrtkappat,w3)
Rho = diag(rhot(t));
Rinv = diag(1./rhot(t));
K = diag(sqrtkappat(t));
W3 = diag(w3(t));
out = delta*kappab/rhob*inv(diag(vol))*K*Rho*C*K*Rinv + W3;
end




