%%% Computes the quasifrequencies for resonator-modulated systems with
%%% a 1D linear lattice inside 2D space

clear all;
close all;
fsz = 16;
lfsz = 14;
alw = 0.75;
lw = 1.2;
N=2;
R=2/(5*N);
epsr=0.2;
epsk=0;
Omega=0.2;
pointM = pi;



    %%% Width of unit cell
    L = 1; 

    %%% Set reciprocal vectors
    b1 = pi/L;

    %%% Set symmetry points in reciprocal space
    pointG = 0;

    %%% Radius of bubble
    R=2/(5*N);
    vol = pi*R^2;

    %%% Material parameters for bubble
    rho_0 = 9000;
    rho_b = 1;
    kappa_0 = 9000;
    kappa_b = 1;

    %%% High contrast parameter \delta
    delta = rho_b/rho_0;
    weight = delta*kappa_b/rho_b/vol;

    %%% Set positions of  bubbles in unit cell
    %%% Dimer
    l = 0.5*R*[1,0];
    c =[];
    for i = 1:N
        c=[c;3*l+((5*l)*(i-1))];
    end
    %%% Define the time modulation
    rhot = @(t) [1./(1+epsr*cos(Omega*t+[0:1/(N-1):1].*pi))];
    % sqrtkappat = @(t) [sqrt(1+epsk*cos(Omega*t)); sqrt(1+epsk*cos(Omega+0.5*pi)); sqrt(1+epsk*cos(Omega+pi))];
    % %w3 = @(t) [3/4*(epsk*cos(Omega*t)/(1+epsk*sin(Omega*t))).^2 + epsk*sin(Omega*t)/(1+epsk*sin(Omega*t));3/4*(epsk*cos(Omega*t+pi)/(1+epsk*sin(Omega*t+pi))).^2 + epsk*sin(Omega*t+pi)/(1+epsk*sin(Omega*t+pi))]; % w3 = -sqrt(k)/2*d/dt
    % w3 = @(t) [epsk*Omega^2*(-4*cos(Omega*t) +epsk*(-5+cos(2*(Omega*t))) ) / (8*(1 + epsk*cos(Omega*t))^2); epsk*Omega^2*(-4*cos(Omega*t+0.5*pi) +epsk*(-5+cos(2*(Omega*t+0.5*pi))) ) / (8*(1 + epsk*cos(Omega*t+0.5*pi))^2); epsk*Omega^2*(-4*cos(Omega*t+pi) +epsk*(-5+cos(2*(Omega*t+pi))) ) / (8*(1 + epsk*cos(Omega*t+pi))^2)];


    sqrtkappat = @(t) [1./sqrt(1+epsk*cos(Omega*t+[0:1/(N-1):1].*pi))];
    w3 = @(t) [epsk*Omega^2*(4*cos(Omega*t+[0:1/(N-1):1]*pi)+epsk*(3+cos(2*(Omega*t+[0:1/(N-1):1]*pi))))./(8*(1+epsk*cos(Omega*t+[0:1/(N-1):1].*pi)).^2)];
    T = 2*pi/Omega;

    %%% Set the discretization or truncation parameters
    N_multipole=3;
    N_a = 100;
    alps = linspace(-pointM,pointM,N_a);

    %%% Precompute values of zeta function (used in computation of lattice sums, just for acceleration)
    d_zeta=makezetadata;
    BandData = zeros(2*N,N_a);
    cnd = zeros(1,N_a);

    k0 = 0.001;
    JHdata = makeJHdata0(k0,R,N_multipole);
    JHijdata = makeJHijexpdata(k0,c,N_multipole);

    alp_ind = 50;
     %%% Set quasi-periodic parameter \alpha
    alp =  alps(alp_ind);
    C = makeC_1D(k0,R,alp,L,d_zeta,JHdata,JHijdata,N,N_multipole);
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    tspan=linspace(0,2);
    odefun = @(t,y) [zeros(2) -M(t); eye(2) zeros(2)]*y;
    y1 = [1 ;0 ;0; 0];
    y2 = [0 ;1 ;0; 0];
    y3= [0 ;0 ;1; 0];
    y4= [0 ;0 ;0; 1];
    [tout,y1out]=ode45(odefun,tspan,y1);
    [tout,y2out]=ode45(odefun,tspan,y2);
    [tout,y3out]=ode45(odefun,tspan,y3);
    [tout,y4out]=ode45(odefun,tspan,y4);
    hold on
    plot(tout,imag(y1out))
    plot(tout,imag(y2out))
    plot(tout,imag(y3out))
    plot(tout,imag(y4out))
    legend(string(1:4))
    %% Plots the solution and saves good images
    
    [w_real_temp, I_sort] = sort(real(BandData)); % sort(real(bands),1);


    w_real = w_real_temp;

    for I_a=1:N_a
        for I_w=1:N
            if w_real(I_w,I_a) > Omega/2*(1-1e-6)
                w_real(I_w,I_a) = w_real(I_w,I_a) - Omega;
            end
        end
    end
    
    for I_a=1:N_a
        for I_w=1:2*N
        w_imag(I_w,I_a) = imag(BandData(I_sort(I_w,I_a),I_a)); % sort(imag(bands),1);
        end
    end
    w_imag = sort(w_imag(1:2*N,:));



%%% Assemble matrix M
function out = Mfunc(t,delta,kappab,rhob,vol,C,rhot,sqrtkappat,w3)
R = diag(rhot(t));
Rinv = diag(1./rhot(t));
K = diag(sqrtkappat(t));
W3 = diag(w3(t));
out = delta*kappab/(vol*rhob)*K*R*C*K*Rinv + W3;
end