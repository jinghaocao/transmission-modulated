%%% This is used to produce k-gap for square lattice
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

% Pick an alpha 
alp = 1/2*(pointM+pointK);


% plot static case
eps = 0; 
epsk = 0;
Omega = 0.2;
modr = @(t,phi) 1./(1+eps*cos(Omega*t+phi*(1)));
rhot = @(t) [modr(t,0+randn(1)); modr(t,2*pi/3+randn(1)); modr(t,4*pi/3+randn(1)); modr(t,0+randn(1)); modr(t,2*pi/3+randn(1)); modr(t,4*pi/3+randn(1))];
modk = @(t,phi) sqrt(1+epsk*cos(Omega*t+phi));
phi_kappa = [0 2*pi/3 4*pi/3 0 4*pi/3 2*pi/3];
sqrtkappat = @(t) [1./sqrt(1+cos(Omega*t+phi_kappa).*epsk)];
w3 = @(t) [Omega^2/4*(1+((epsk.^2-1)./(1+epsk.*cos(Omega*t+phi_kappa)).^2))];
T = 2*pi/Omega;

% Compute the bloch modes at that alpha
C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
ones_t = @(t) ones(N,1);
zeros_t = @(t) zeros(N,1);
M0= @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,ones_t,ones_t,zeros_t);
W_static = make_bloch_mode(M,M0,N,T);
figure
    title('static case')

for N_mode = 1:N;
    subplot(2,N/2,N_mode) 
    hold on
    t = linspace(0,2*pi);
    for i = -2:2
        for j = -2:2
        % Quasiperiodicity 
        color_array =real(W_static(:,N_mode)*exp(1i*dot(alp,[i*L1(1)+j*L2(1),i*L1(2)+j*L2(2)])));
        for n = 1:N
            cc = color_array(n);
            if cc > 0
                scatter(c(n,1)+i*L1(1)+j*L2(1),c(n,2)+i*L1(2)+j*L2(2),50,[1-cc 1-cc 1-cc],'filled')
            else
                scatter(c(n,1)+i*L1(1)+j*L2(1),c(n,2)+i*L1(2)+j*L2(2),50,[0 1+cc 0],'filled')
            end
        end
        end
    end
    daspect([1 1 1])
end

% introduce randomness for time modulation

%%% Define the time modulation
eps = 0.2; 
epsk = 0;
Omega = 0.2;
W = zeros(2*N,2*N);
for r = 1:100
    modr = @(t,phi) 1./(1+eps*cos(Omega*t+phi*(1+randn(1))));
    rhot = @(t) [modr(t,0); modr(t,2*pi/3); modr(t,4*pi/3); modr(t,0); modr(t,2*pi/3); modr(t,4*pi/3)];
    modk = @(t,phi) sqrt(1+epsk*cos(Omega*t+phi));
    phi_kappa = [0 2*pi/3 4*pi/3 0 4*pi/3 2*pi/3];
    sqrtkappat = @(t) [1./sqrt(1+cos(Omega*t+phi_kappa).*epsk)];
    w3 = @(t) [Omega^2/4*(1+((epsk.^2-1)./(1+epsk.*cos(Omega*t+phi_kappa)).^2))];
    T = 2*pi/Omega;

    % Compute the bloch modes at that alpha
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    ones_t = @(t) ones(N,1);
    zeros_t = @(t) zeros(N,1);
    M0= @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,ones_t,ones_t,zeros_t);
    W = W+make_bloch_mode(M,M0,N,T);
end
W = W/100;
% Plot the mode in 2D
figure
for N_mode = 1:N;
    subplot(2,N/2,N_mode) 
    hold on
    t = linspace(0,2*pi);
    for i = -2:2
        for j = -2:2
        % Quasiperiodicity 
        color_array =real(W(:,N_mode)*exp(1i*dot(alp,[i*L1(1)+j*L2(1),i*L1(2)+j*L2(2)])));
        for n = 1:N
            cc = color_array(n);
            if cc > 0
                scatter(c(n,1)+i*L1(1)+j*L2(1),c(n,2)+i*L1(2)+j*L2(2),50,[1-cc 1-cc 1-cc],'filled')
            else
                scatter(c(n,1)+i*L1(1)+j*L2(1),c(n,2)+i*L1(2)+j*L2(2),50,[0 1+cc 0],'filled')
            end
        end
        end
    end
    daspect([1 1 1])
end
    title('randomlized with 100 repititions,epsr=0.2,epsk=0')


%%% Assemble matrix M
function out = Mfunc(t,delta,kappab,rhob,vol,C,rhot,sqrtkappat,w3)
    Rho = diag(rhot(t));
    Rinv = diag(1./rhot(t));
    K = diag(sqrtkappat(t));
    W3 = diag(w3(t));
    out = delta*kappab/rhob*inv(diag(vol))*K*Rho*C*K*Rinv + W3;
end