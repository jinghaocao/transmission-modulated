%%% This is used to produce k-gap for honeycomb lattice
%%% when modulating kappa only. Figure 7

%%% lattice vectors ( a1 = [L1x;0], a2 = [L2x;L2y], D: distance between nearest neighbor bubbles ) 
D = 1; 
L1x = D*sqrt(3);
L2x = D*sqrt(3)/2;
L2y = D*3/2;
L1 = [L1x,0];
L2 = [L2x,L2y];

%%% Set reciprocal vectors
b1 = 4*pi/(3*D)*[sqrt(3)/2, -1/2];
b2 = 4*pi/(3*D)*[0, 1];

%%% Set symmetry points in reciprocal space
pointG = [0,0];
pointM = 1/2*b1;
pointK = 2/3*b1 + 1/3*b2;

%%% Radius of bubble
R=D*0.1; %R = D*0.1;
vol = pi*R^2;

%%% Material parameters for bubble
rho_0 = 9000;
rho_b = 1;
kappa_0 = 9000;
kappa_b = 1;

%%% High contrast parameter \delta
delta = rho_b/rho_0;
weight = delta*kappa_b/rho_b/vol;

%%% Set positions of two bubbles in unit cell

%%% Standard honeycomb
N = 2; % Number of resonators in unit cell
c1 = 1/3*(L1+L2);
c2 = 2/3*(L1+L2);
c = [c1;c2];

%%% Hexagon with trimer vertices
N = 6; % Number of resonators in unit cell
c1 = 1/3*(L1+L2);
c2 = 2/3*(L1+L2);
eps = 3*R; %eps = 3*R;
l = @(theta) eps*[cos(theta+pi/6),sin(theta+pi/6)];
c = [c1+l(0);c1+l(2*pi/3);c1+l(4*pi/3);c2+l(pi/3);c2+l(3*pi/3);c2+l(5*pi/3)];

%%% Extended or shrunk hexagons
% N = 6; % Number of resonators in unit cell
% c0 = 0.5*(L1+L2);
% theta0 = atan2(c0(2),c0(1))+pi/2;
% eps = 1.1; 
% l = @(theta) eps*norm(L1-L2)/3*[cos(theta+theta0),sin(theta+theta0)];
% c = [c0+l(0);c0+l(pi/3);c0+l(2*pi/3);c0+l(3*pi/3);c0+l(4*pi/3);c0+l(5*pi/3)];

%%% Plot the geometry
figure, hold on
t = linspace(0,2*pi);
for i = -2:2
    for j = -2:2
        for n = 1:N
            plot(c(n,1)+R*cos(t)+i*L1(1)+j*L2(1), c(n,2)+R*sin(t)+i*L1(2)+j*L2(2),'k')  
            text(c(n,1),c(n,2),num2str(n))
        end
    end
end
daspect([1 1 1])
hold off
close
%%% Define the time modulation
eps = 0; 
epsk = 0.2;
Omega = 0.2;

modr = @(t,phi) 1./(1+eps*cos(Omega*t+phi));
rhot = @(t) [modr(t,0); modr(t,2*pi/3); modr(t,4*pi/3); modr(t,0); modr(t,2*pi/3); modr(t,4*pi/3)];
modk = @(t,phi) sqrt(1+epsk*cos(Omega*t+phi));
phi_kappa = [0 2*pi/3 4*pi/3 0 4*pi/3 2*pi/3];
sqrtkappat = @(t) [1./sqrt(1+cos(Omega*t+phi_kappa).*epsk)];
w3 = @(t) [Omega^2/4*(1+((epsk.^2-1)./(1+epsk.*cos(Omega*t+phi_kappa)).^2))];
T = 2*pi/Omega;
%%% Set the discretization or truncation parameters
N_multipole=3;
N_lattice=30;
N_KM = 100;
N_MG = floor(sqrt(3)*N_KM);
N_GK = 2*N_KM;
%%% Precompute values of zeta function (used in computation of lattice sums, just for acceleration)
d_zeta=makezetadata;
BandDataKM = zeros(2*N,N_KM);
BandDataMG = zeros(2*N,N_MG);
BandDataGK = zeros(2*N,N_GK);
k0 = 0.001;
JHdata = makeJHdata0original(k0,R,N_multipole);
JHijdata = makeJHijexpdata(k0,c,N_multipole);
    
for alp_ind = 1:N_MG
    %%% Set quasi-periodic parameter \alpha
    alp =  pointM + (alp_ind-1)/N_MG * (pointG-pointM);
    %%% Compute the capacitance matrix C
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    BandDataMG(:, alp_ind) = hill_exp(T,M,N);
end

parfor alp_ind = 1:N_GK    
    %%% Set quasi-periodic parameter \alpha
    alp =  pointG + (alp_ind-1)/N_GK * (pointK-pointG);
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    BandDataGK(:, alp_ind) = hill_exp(T,M,N);
end

parfor alp_ind = 1:N_KM
    %%% Set quasi-periodic parameter \alpha
    alp =  pointK + (alp_ind-1)/N_KM * (pointM-pointK);

    %%% Compute the operator A
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    BandDataKM(:, alp_ind) = hill_exp(T,M,N);
end    

lengthMG = norm(pointG-pointM);
lengthGK = norm(pointK-pointG);
lengthKM = norm(pointM-pointK);
lengthMGK = lengthMG + lengthGK + lengthKM;
alps = [(0:N_MG-1)*lengthMG/N_MG, (0:N_GK-1)*lengthGK/N_GK + lengthMG, (0:N_KM-1)*lengthKM/N_KM + lengthMG + lengthGK];
a_plot = [(0:N_MG-1)*lengthMG/N_MG, (0:N_GK-1)*lengthGK/N_GK + lengthMG, (0:N_KM-1)*lengthKM/N_KM + lengthMG + lengthGK];
w_plot = [BandDataMG, BandDataGK, BandDataKM];
N_a = length(a_plot);
% Plots the solution and saves good images
fsz = 16;
lfsz = 14;
alw = 0.75;
lw = 1.2;
w_real = sort(real(w_plot),1);
w_imag = sort(imag(w_plot),1);
for I_a=1:N_a
    for I_w=N+1:2*N
        if w_real(I_w,I_a) > Omega/2*(1-1e-3)
            w_real(I_w,I_a) = w_real(I_w,I_a) - Omega;
        end
    end
end
color_array = {'k','b','r','g','y',[.5 .6 .7]};
figure
hold on
for I_w = N+1:2*N
plot(a_plot,w_real(I_w,:),'.b','MarkerSize',2.5)
end
hold on
for I_w = N+1:2*N
plot(a_plot,w_imag(I_w,:),'.r')
end

x = [ 0*pi, lengthMG ,  lengthMG+lengthGK,      lengthMGK];
y= ['     M  ';'$\Gamma$';'    K   ';'   M    '];
axis([0,lengthMGK,0,Omega/2])
set(gca,'XTick',x); % Change x-axis ticks
set(gca,'XTickLabel',y); 
set(gca,'fontsize', 15);
set(gca,'xgrid','on')
%legend('Real part','Imaginary part','position',[0.73,0.63,0.2,0.1],'interpreter','latex');
xlabel('Quasi-periodicity $\alpha$','interpreter','latex');
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
ax.Position = [left bottom ax_width ax_height*0.99];
show_variable = strcat('epsr =',num2str(eps),'epsk=',num2str(epsk));
% Uncomment to plot k-gaps
line([0,0],[0,0.0876074],'Color','g');
line([0.266339,0.266339],[0,0.0888022],'Color','g');
line([0.399509,0.399509],[0,0.060091],'Color','g');
line([0.859549,0.859549],[0,0.0529812],'Color','g');
line([3.33987,3.33987],[0,0.0525949],'Color','g');
line([3.69054,3.69054],[0,0.0585829],'Color','g');
line([3.66635,3.66635],[0,0.0879411],'Color','g');
line([4.36769,4.36769],[0,0.0843817],'Color','g');
line([4.66999,4.66999],[0,0.0852427],'Color','g');
line([5.7099,5.7099],[0,0.0876074],'Color','g');

% Choose the right legend
lgd =legend('Re(\omega)','','','','','','','Im(\omega)','','','','','k-gap')
lgd =legend('Re(\omega)','','','','','','','Im(\omega)','','','')
lgd.Location = 'east';
lgd.FontSize = 8;

%%% Assemble matrix M
function out = Mfunc(t,delta,kappab,rhob,vol,C,rhot,sqrtkappat,w3)
R = diag(rhot(t));
Rinv = diag(1./rhot(t));
K = diag(sqrtkappat(t));
W3 = diag(w3(t));
out = delta*kappab/(vol*rhob)*K*R*C*K*Rinv + W3;
end
