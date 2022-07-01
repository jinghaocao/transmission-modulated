
 %%% lattice vectors ( a1 = [L1x;0], a2 = [L2x;L2y], D: distance between nearest neighbor bubbles ) 
 Omega = 0.2  
 D = 1; 
    L1x = D*sqrt(3);
    L2x = D*sqrt(3)/2;
    L2y = D*3/2;
    L1 = [L1x,0];
    L2 = [L2x,L2y];

%     %%% Set reciprocal vectors
%     b1 = 4*pi/(3*D)*[sqrt(3)/2, -1/2];
%     b2 = 4*pi/(3*D)*[0, 1];

%     %%% Set symmetry points in reciprocal space
%     pointG = [0,0];
%     pointM = angle*1/2*b1;
%     pointK = angle*(2/3*b1 + 1/3*b2);

    %%% Radius of bubble
    R=0.1; %R = D*0.1;
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



D =1;
%%% Set reciprocal vectors
x = sqrt(3);
b1 = 4*pi/(3*D)*[x/2, -1/2];
b2 = 4*pi/(3*D)*[0, 1];

%%% Set symmetry points in reciprocal space
pointG = [0,0];
pointM = 1/2*b1;
pointK = -2/3*b1 - 1/3*b2;

alpha_wanted = pointK;

  epsk = 0;
    modr = @(t,phi) 1./(1+epsr*cos(Omega*t+phi));
    rhot = @(t) [modr(t,0); modr(t,2*pi/3); modr(t,4*pi/3); modr(t,0); modr(t,2*pi/3); modr(t,4*pi/3)];
    phi = [0 2*pi/3 4*pi/3 0 2*pi/3 4*pi/3];
    modk = @(t,phi) sqrt(1+epsk*cos(Omega*t+phi));
    w3k = @(t,phi) epsk*Omega^2*(-4*cos(Omega*t+phi) +epsk*(-5+cos(2*(Omega*t+phi))) ) / (8*(1 + epsk*cos(Omega*t+phi))^2);
    %w3k = @(t,phi) -3/4*(epsk*cos(Omega*t+phi)/(1+epsk*sin(Omega*t+phi))).^2 + epsk*sin(Omega*t+phi)/(1+epsk*sin(Omega*t+phi));
    sqrtkappat = @(t) [modk(t,0); modk(t,2*pi/3); modk(t,4*pi/3); modk(t,0); modk(t,2*pi/3); modk(t,4*pi/3)];
    w3 = @(t)  [w3k(t,0); w3k(t,2*pi/3); w3k(t,4*pi/3); w3k(t,0); w3k(t,2*pi/3); w3k(t,4*pi/3)];
    TT = 2*pi/Omega;
    
N_multipole=4;
N_a = 1000;
% alps = linspace(-pointM,pointM,N_a);
%%% Precompute values of zeta function (used in computation of lattice sums, just for acceleration)
d_zeta=makezetadata;
% BandData = zeros(2*N,N_a);
% cnd = zeros(1,N_a);

k0 = 0.00001;
[Jdata,Hdata] = makeJHdata0(k0,R,N_multipole);
JHijdata = makeJHijexpdata(k0,c,N_multipole);
p_plus = makepvalue(alpha_wanted,delta,kappa_b,vol,Omega,phi)
p_minus = makepvalue(-alpha_wanted,delta,kappa_b,vol,Omega,phi)

function p = makepvalue(alpha_wanted,delta,kappa_b,vol,Omega,phi)

    TT = 2*pi/Omega;
%%% lattice vectors ( a1 = [L1x;0], a2 = [L2x;L2y], D: distance between nearest neighbor bubbles ) 
    D = 1; 
    L1x = D*sqrt(3);
    L2x = D*sqrt(3)/2;
    L2y = D*3/2;
    L1 = [L1x,0];
    L2 = [L2x,L2y];

%     %%% Set reciprocal vectors
%     b1 = 4*pi/(3*D)*[sqrt(3)/2, -1/2];
%     b2 = 4*pi/(3*D)*[0, 1];

%     %%% Set symmetry points in reciprocal space
%     pointG = [0,0];
%     pointM = angle*1/2*b1;
%     pointK = angle*(2/3*b1 + 1/3*b2);

    %%% Radius of bubble
    R=0.1; %R = D*0.1;
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

%     Plot the geometry
%     figure, hold on
%     t = linspace(0,2*pi);
%     for i = -2:2
%         for j = -2:2
%             for n = 1:N
%                 plot(c(n,1)+R*cos(t)+i*L1(1)+j*L2(1), c(n,2)+R*sin(t)+i*L1(2)+j*L2(2),'k')  
%                 text(c(n,1),c(n,2),num2str(n))
%             end
%         end
%     end
%     daspect([1 1 1])
%     hold off

    %%% Define the time modulation
    epsk = 0;
    epsr = 0;
    modr = @(t,phi) 1./(1+epsr*cos(Omega*t+phi));
    rhot = @(t) [modr(t,0); modr(t,2*pi/3); modr(t,4*pi/3); modr(t,0); modr(t,2*pi/3); modr(t,4*pi/3)];
    modk = @(t,phi) sqrt(1+epsk*cos(Omega*t+phi));
    w3k = @(t,phi) epsk*Omega^2*(-4*cos(Omega*t+phi) +epsk*(-5+cos(2*(Omega*t+phi))) ) / (8*(1 + epsk*cos(Omega*t+phi))^2);
    %w3k = @(t,phi) -3/4*(epsk*cos(Omega*t+phi)/(1+epsk*sin(Omega*t+phi))).^2 + epsk*sin(Omega*t+phi)/(1+epsk*sin(Omega*t+phi));
    sqrtkappat = @(t) [modk(t,0); modk(t,2*pi/3); modk(t,4*pi/3); modk(t,0); modk(t,2*pi/3); modk(t,4*pi/3)];
    w3 = @(t)  [w3k(t,0); w3k(t,2*pi/3); w3k(t,4*pi/3); w3k(t,0); w3k(t,2*pi/3); w3k(t,4*pi/3)];
    T = 2*pi/Omega;

    %%% Set the discretization or truncation parameters
    N_multipole=3;
    N_lattice=30;



    %%% Precompute values of zeta function (used in computation of lattice sums, just for acceleration)
    d_zeta=makezetadata;

    k0 = 0.001;
    JHdata = makeJHdata0original(k0,R,N_multipole);
    JHijdata = makeJHijexpdata(k0,c,N_multipole);
    
    
    alp =  alpha_wanted;
    C = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multipole,N_lattice);
    M_t = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    M=M_t(0);
    A= [zeros(N) eye(N); -M_t(0) zeros(N)];
    [T,D]=eig(A,'vector');
    F0=zeros(N*2);
    m = round(imag(D)/Omega)
    F0 = real(D)+1i*(imag(D)-(round(imag(D)/Omega))*Omega);
    f0 = real(F0/(1i));
    [f0,ind]= sort(f0,'descend');
    m = m(ind);
    T= T(:,ind);
    M1 = zeros(N);
    M2 = zeros(N);
    for l = 1 : N
        for j = 1: N
            M1(l,j)=-0.5*M(l,j)*(exp(1i*phi(l))-exp(1i*phi(j)));
            M2(l,j)=-0.5*M(l,j)*(exp(-1i*phi(l))-exp(-1i*phi(j)));
        end
    end
if m(1)-m(2)==1
    A1_1 = inv(T)*[zeros(N) zeros(N); M1 zeros(N)]*T;
    A1_n1 = inv(T)*[zeros(N) zeros(N); M2 zeros(N)]*T;
end
if m(1)-m(2)==-1
    A1_n1 = inv(T)*[zeros(N) zeros(N); M1 zeros(N)]*T;
    A1_1 = inv(T)*[zeros(N) zeros(N); M2 zeros(N)]*T; 
end
    p = imag(sqrt(A1_1(1,2)*A1_n1(2,1)));
end
    
function out = Mfunc(t,delta,kappab,rhob,vol,C,rhot,sqrtkappat,w3)
Rho = diag(rhot(t));
Rinv = diag(1./rhot(t));
K = diag(sqrtkappat(t));
W3 = diag(w3(t));
out = delta*kappab/rhob*inv(diag(vol))*K*Rho*C*K*Rinv + W3;
end