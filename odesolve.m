%%% Computes the quasifrequencies for resonator-modulated systems with
%%% a 1D linear lattice inside 2D space
close all
clear all

fsz = 16;
lfsz = 14;
alw = 0.75;
lw = 1.2;
N=2;
R=2/(5*N);
epsr=0.1;
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
    rhot = @(t) [1./(1+epsr*cos(Omega*t+[0:1/(N-1):1].*0.3*pi))];
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

     %%% Set quasi-periodic parameter \alpha
    alp =  -1;
    C = makeC_1D(k0,R,alp,L,d_zeta,JHdata,JHijdata,N,N_multipole);
    M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
    
    W = zeros(2*N,2*N);
    Z = zeros(N,N);
    I = eye(N,N);
    MM = @(t,y) [Z, I; -M(t), Z]*y;
    for j = 1:N
        wI0 = zeros(2*N,1); wI0(j) = 1;
        wII0 = zeros(2*N,1); wII0(N+j) = 1;
        [~, wI] = ode45(MM, [0,T], wI0);
        [~, wII] = ode45(MM, [0,T], wII0);
        W(:,j) = wI(end,:).';
        W(:,N+j) = wII(end,:).';
    end
    
    for j = 1:N
        wI0 = zeros(2*N,1); wI0(j) = 1;
        wII0 = zeros(2*N,1); wII0(N+j) = 1;
        [~, wI] = ode45(MM, [0,-10*T], wI0);
        [~, wII] = ode45(MM, [0,-10*T], wII0);
        Wminus(:,j) = wI(end,:).';
        Wminus(:,N+j) = wII(end,:).';
    end
    %% Floquet exponent
    FE=real(log(eig(W))/(1i*T));
    FE
    W
    Wminus
 
    %%  test claim 1
    P= [-eye(N),zeros(N);zeros(N),eye(N)];
    W-P*Wminus*P
    
    %%  test claim 2
    W*Wminus
    
    
    %%  test Theorem
    W*P*W

    %% Plots the solution and saves good images
    
    parfor alp_ind = 1:N_a
        %%% Set quasi-periodic parameter \alpha
        alp2 =  alps(alp_ind);

        %%% Compute the operator A
        C = makeC_1D(k0,R,alp2,L,d_zeta,JHdata,JHijdata,N,N_multipole);
        M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
        [BandData(:, alp_ind), cnd(alp_ind)] = hill_exp(T,M,N);
    end
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
    hold on
    for I_w = 1:N*2
        plot(alps,w_real(I_w,:),'b-','LineWidth',1)
        plot(-alps,w_real(I_w,:),'r--','LineWidth',2)
    end
    
    legend('alpha','-alpha');
x = [ -pointM, 0,  pointM];
y= ['$-\pi/L$';'$\Gamma$';' $\pi/L$'];
xlim([-pointM,pointM])
set(gca,'XTick',x); % Change x-axis ticks
set(gca,'XTickLabel',y); 

set(gca,'xgrid','on')

%legend('Real part','Imaginary part','position',[0.74,0.62,0.2,0.1],'interpreter','latex','FontSize',lfsz);
xlabel('Quasiperiodicity $\alpha$','interpreter','latex');
ylabel('Frequency $\omega$','interpreter','latex');
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
%print('SLresmod_b','-depsc');
print('SL_c_def','-depsc');

%%% Assemble matrix M
function out = Mfunc(t,delta,kappab,rhob,vol,C,rhot,sqrtkappat,w3)
R = diag(rhot(t));
Rinv = diag(1./rhot(t));
K = diag(sqrtkappat(t));
W3 = diag(w3(t));
out = delta*kappab/(vol*rhob)*K*R*C*K*Rinv + W3;
end