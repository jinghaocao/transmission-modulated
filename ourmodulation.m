%%% Computes the quasifrequencies for resonator-modulated systems with
%%% a 1D linear lattice inside 2D space
clear all;
close all;



fsz = 16;
lfsz = 14;
alw = 0.75;
lw = 1.2;

%%% Width of unit cell
L = 1; 

%%% Set reciprocal vectors
b1 = pi/L;

%%% Set symmetry points in reciprocal space
pointG = 0;
pointM = b1;

%%% Radius of bubble
R=L*[0.1; 0.1;0.1];
vol = pi*R.^2;

%%% Material parameters for bubble
rho_0 = 9000;
rho_b = 1;
kappa_0 = 9000;
kappa_b = 1;

%%% High contrast parameter \delta
delta = rho_b/rho_0;
% weight = delta*kappa_b/rho_b/vol;

%%% Set positions of two bubbles in unit cell
%%% Dimer

N = size(R,1); % Number of resonators in unit cell
l = 0.95/(5*N)*[1,0];
c =[];
for i = 1:N
    c=[c;3*l+((5*l)*(i-1))];
end

%%% Plot the geometry
% figure, hold on
% t = linspace(0,2*pi);
% for i = -2:2
%     for n = 1:N
%         plot(c(n,1)+R(n)*cos(t)+i*L, c(n,2)+R(n)*sin(t),'k')  
%     end
% end
% daspect([1 1 1])
% hold off

%%% Define the time modulation
Omega = 0.3; %Omega = 0.3; gives us band gap
alpha_wanted = -0.32;

phi = [0 pi/2 pi];
epsr = 0;
rhot = @(t) [1./(1+(cos(Omega*t+phi).*epsr))];
epsk= 0;
sqrtkappat = @(t) [1./sqrt(1+cos(Omega*t+phi)*epsk)];
w3 = @(t) [Omega^2*(4*cos(Omega*t+[0:1/(N-1):1]*0.5*pi)*epsk+(3+cos(2*(Omega*t+[0:1/(N-1):1]*0.5*pi)))*epsk)./(8*(1+cos(Omega*t+[0:1/(N-1):1].*0.5*pi)*epsk).^2)];
T = 2*pi/Omega;
N_a = 150;

%% get static band data
[w_real,w_imag,alps]= makeBandData(pointM,c,Omega,delta,kappa_b,rho_b,R,L,vol,rhot,sqrtkappat,w3,N_a);

%% Plots the solution and saves good images
fsz = 16;
lfsz = 14;
alw = 0.75;
lw = 1.2;

figure
hold on
for I_a=1:N_a
    for I_w=1:N*2
        if w_real(I_w,I_a) > Omega/2*(1-1e-6)
            w_real(I_w,I_a) = w_real(I_w,I_a) - Omega;
        end
    end
end

for I_w = N+1:N*2
    h1 = plot(alps,w_real(I_w,:),'.b','LineWidth',lw)
  
    h2 = plot(alps,w_imag(I_w,:),'.m','LineWidth',lw)
% 
%     if I_w< N*2 & max(w_real(I_w,1:N_a/2)) < min(w_real(I_w+1,1:N_a/2))
%         h3 = plot(alps(1:N_a/2),ones(N_a/2,1)*max(w_real(I_w,1:N_a/2)),'r','LineWidth',lw);
%         plot(alps(1:N_a/2),ones(N_a/2,1)*min(w_real(I_w+1,1:N_a/2)),'r','LineWidth',lw);
%     end
%     if I_w< N*2 & max(w_real(I_w,N_a/2+1:end)) < min(w_real(I_w+1,N_a/2+1:end))
%         plot(alps(N_a/2+1:end),ones(N_a/2,1)*max(w_real(I_w,N_a/2+1:end)),'r','LineWidth',lw);
%         plot(alps(N_a/2+1:end),ones(N_a/2,1)*min(w_real(I_w+1,N_a/2+1:end)),'r','LineWidth',lw);
%     end
    

%     plot(-alps,w_real(I_w,:),'.r','LineWidth',lw)
%     plot(alps,w_imag(I_w,:),'-.g','LineWidth',lw)
%     plot(-alps,w_imag(I_w,:),'.r','LineWidth',lw)

end
legend([h1,h2],'Real part','Imaginary part','Location','best'	);

% legend([h1,h2,h3],'Real part','Imaginary part','Bandgap','Location','best'	);
x = [ -pointM, 0,  pointM];
y= ['$-\pi/L$';'$\Gamma$';' $\pi/L$'];
ylim([-0.02,0.15])
xlim([-pointM,pointM])
set(gca,'XTick',x); % Change x-axis ticks
set(gca,'XTickLabel',y); 

set(gca,'xgrid','on')
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
print('epsilon0','-depsc');

% Precompute values of zeta function (used in computation of lattice sums, just for acceleration)
% d_zeta=makezetadata;
% % BandData = zeros(2*N,N_a);
% % cnd = zeros(1,N_a);
% 
% k0 = 0.00001;
% N_multipole=4;
% [Jdata,Hdata] = makeJHdata0(k0,R,N_multipole);
% JHijdata = makeJHijexpdata(k0,c,N_multipole);
% alpha_wanted = compute_degeneracy(c,R,L,Omega,kappa_b,rho_b,vol,delta);
% alpha_wanted
% %Plot perturbation at degenercy point
% epsilon = [0:0.001:0.5];
% for epsilon_index = 1:length(epsilon)
%     epsr = epsilon(epsilon_index);
%     rhot = @(t) [1./(1+(cos(Omega*t+phi)*epsr))];
%     C = makeC_1D(k0,R,alpha_wanted,L,d_zeta,Jdata,Hdata,JHijdata,N,N_multipole);
%     M_1 = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
%     [BandData2(:, epsilon_index), cnd] = hill_exp(T,M_1,N);
%     BandDatareal(:, epsilon_index)=sort(real(BandData2(:, epsilon_index)));
%     
%     C_2 = makeC_1D(k0,R,-alpha_wanted,L,d_zeta,Jdata,Hdata,JHijdata,N,N_multipole);
%     M_2 = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C_2,rhot,sqrtkappat,w3);
%     [BandData3(:, epsilon_index), cnd] = hill_exp(T,M_2,N);
%     BandDatareal2(:, epsilon_index)=sort(real(BandData3(:, epsilon_index)));
% end
% 
% figure
% poly = zeros(6,5);
% for i = 6 :-1: 5
%     hold on
%     plot(epsilon,BandDatareal(i,:),'r','LineWidth',1);
%     hold on
%     plot(epsilon,BandDatareal2(i,:),'b','LineWidth',1.2);
% end
% 
% ax = gca;
% set(gca,'TickLabelInterpreter','latex')
% set(gca, 'FontSize', fsz, 'LineWidth', alw);
% set(gca,'TickLabelInterpreter','latex')
% 
% hold on
%  plot(epsilon,BandDatareal(4,:),'r','LineWidth',1);
%     plot(epsilon,BandDatareal2(4,:),'b','LineWidth',1.2);
%      xlabel('Modulation parameter $\varepsilon$','interpreter','latex')
%     ylabel('Non-degenerate frequencies','interpreter','latex')
%     set(gca, 'FontSize', fsz, 'LineWidth', alw);
% set(gca,'TickLabelInterpreter','latex')
% legend(strcat('$\alpha =$ ' , num2str(alpha_wanted)),strcat('$\alpha =$ ' , num2str(-alpha_wanted)),'interpreter','latex','Location','southwest')
% 
% j =1;
% for i = 6 :-1: 4
%     poly(j,:) = polyfit(epsilon,BandDatareal(i,:),4);
%     poly(j+1,:) = polyfit(epsilon,BandDatareal2(i,:),4);
%     j = j + 2;
% end
