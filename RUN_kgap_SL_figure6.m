%%% This is used to produce k-gap for square lattice
%%% when modulating kappa only. Figure 7


%%% lattice vectors ( a1 = [L1x;0], a2 = [L2x;L2y], D: distance between nearest neighbor bubbles ) 
D = 1; 
L1x = D;
L2x = 0;
L2y = D;
L1 = [L1x,0];
L2 = [L2x,L2y];
%%% Set reciprocal vectors
b1 = 2*pi/L1x*[1, 0];
b2 = 2*pi/L2y*[0, 1];
%%% Set symmetry points in reciprocal space
pointG = [0,0];
pointX = 1/2*b1;
pointM = 1/2*(b1 + b2);
%%% Radius of bubble
R=D*0.1;
vol = pi*R^2;
%%% Material parameters for bubble
rho_0 = 9000;
rho_b = 1;
kappa_0 = 9000;
kappa_b = 1;
%%% High contrast parameter \delta
delta = rho_b/rho_0;
weight = delta*kappa_b/rho_b/vol;
%%% Set positions of 3 bubbles in unit cell
%%% Dimer
N = 3; % Number of resonators in unit cell
c0 = 1/2*(L1+L2);
eps = 3*R;
l = @(theta) eps*[cos(theta+pi/6),sin(theta+pi/6)];
c = [c0+l(0);c0+l(2*pi/3);c0+l(4*pi/3)];
% Plot the geometry
figure, hold on
t = linspace(0,2*pi);
for i = -2:2
    for j = -2:2
        for n = 1:N
            plot(c(n,1)+R*cos(t)+i*L1(1)+j*L2(1), c(n,2)+R*sin(t)+i*L1(2)+j*L2(2),'k')  
        end
    end
end
daspect([1 1 1])
hold off
close
%%% Define the time modulation
epsr = 0;
Omega = 0.2; 
phii = [0 2/3 4/3]*pi;
rhot = @(t) [ 1./(1+epsr*cos(Omega*t+pi))];
epsk = 0.1;
sqrtkappat = @(t) [1./sqrt(1+cos(Omega*t+phii).*epsk)];
w3 = @(t) [Omega^2/4*(1+((epsk.^2-1)./(1+epsk.*cos(Omega*t+phii)).^2))];
T = 2*pi/Omega;
%%% Set the discretization or truncation parameters
N_multipole=3;
N_lattice=10;
N_XG = 100;
N_GM = floor(sqrt(2)*N_XG);
N_MX = 2*N_XG;
%%% Precompute values of zeta function (used in computation of lattice sums, just for acceleration)
d_zeta=makezetadata;
BandDataXG = zeros(2*N,N_XG);
BandDataGM = zeros(2*N,N_GM);
BandDataMX = zeros(2*N,N_MX);
cndXG = zeros(1,N_XG);
cndGM = zeros(1,N_GM);
cndMX = zeros(1,N_MX);
k0 = 0.001;
JHdata = makeJHdata0original(k0,R,N_multipole);
JHijdata = makeJHijexpdata(k0,c,N_multipole);
for alp_ind = 1:N_XG
    %%% Set quasi-periodic parameter \alpha
    alp =  pointX + (alp_ind-1)/N_XG * (pointG-pointX);
    %%% Compute the operator A
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    [BandDataXG(:, alp_ind), cndXG(alp_ind)] = hill_exp(T,M,N);
end  

for alp_ind = 1:N_GM
    %%% Set quasi-periodic parameter \alpha
    alp =  pointG + (alp_ind-1)/N_GM * (pointM-pointG);
    %%% Compute the capacitance matrix C
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    [BandDataGM(:, alp_ind), cndGM(alp_ind)] = hill_exp(T,M,N);
end

parfor alp_ind = 1:N_MX    
    %%% Set quasi-periodic parameter \alpha
    alp =  pointM + (alp_ind-1)/N_MX * (pointX-pointM);
    %%% Compute the operator A
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    [BandDataMX(:, alp_ind), cndMX(alp_ind)] = hill_exp(T,M,N);
end  

lengthXG = norm(pointG-pointX);
lengthGM = norm(pointM-pointG);
lengthMX = norm(pointX-pointM);
lengthXGM = lengthXG + lengthGM + lengthMX;

%% Plots the solution and saves good images
% close all
fsz = 16;
lfsz = 14;
alw = 0.75;
lw = 1.2;
bands = [BandDataXG, BandDataGM, BandDataMX];
alps = [(0:N_XG-1)*lengthXG/N_XG, (0:N_GM-1)*lengthGM/N_GM + lengthXG, (0:N_MX-1)*lengthMX/N_MX + lengthXG + lengthGM];
cdn = [cndXG, cndGM, cndMX];
N_a = length(alps);

w_real = real(bands);
w_imag = imag(bands);

figure
hold on
colorarry = {'b','r','k','m','y',[.8 .2 .6]};

for I_w = 1:N
    plot(alps,w_real(I_w,:),colorarry{1})
    scatter(alps,w_imag(I_w,:),4,colorarry{2},'filled')
end
% Uncommnet to show k-gap when eps_kappa = 0.1
line([0.69115,0.69115],[-Omega/2,Omega/2],'Color','g');
line([5.18973,5.18973],[-Omega/2,Omega/2],'Color','g');
line([5.56785,5.56785],[-Omega/2,Omega/2],'Color','g');
line([9.95638,9.95638],[-Omega/2,Omega/2],'Color','g');

lgd =legend('Re(\omega)','Im(\omega)','','','','','k-gap')
lgd =legend('Re(\omega)','Im(\omega)','','','','')

lgd.FontSize=11;
lgd.Location = 'east'
x = [ 0*pi, lengthXG ,  lengthXG+lengthGM,      lengthXGM];
y= ['     X  ';'$\Gamma$';'    M   ';'   X    '];
ylim([-max(max(w_imag))*1.2,Omega/2])
xlim([0,lengthXGM])
set(gca,'XTick',x); % Change x-axis ticks
set(gca,'XTickLabel',y); 

set(gca,'xgrid','on')

%legend('Real part','Imaginary part','position',[0.74,0.62,0.2,0.1],'interpreter','latex','FontSize',lfsz);
xlabel('Quasiperiodicity $\alpha$','interpreter','latex');
ylabel('Quasifrequency $\omega$','interpreter','latex');
ax = gca;

set(gcf, 'Position', [0, 0, 600, 500])
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(gca,'TickLabelInterpreter','latex')
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left*1.15 bottom ax_width*0.99 ax_height];

%%% Assemble matrix M
function out = Mfunc(t,delta,kappab,rhob,vol,C,rhot,sqrtkappat,w3)
    Rho = diag(rhot(t));
    Rinv = diag(1./rhot(t));
    K = diag(sqrtkappat(t));
    W3 = diag(w3(t));
    out = delta*kappab/rhob*inv(diag(vol))*K*Rho*C*K*Rinv + W3;
end
