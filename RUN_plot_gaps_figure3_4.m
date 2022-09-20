% This is used to produce Figures 3 and 4
Na = 200;
% Change the parameters to suit the purposes
epsr = 0.1;
epsk = 0;

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
%     sort to gain correct band structures, only uncomment in eps_k>0 and
%     eps_rho =0
%     if alp(i)>-1.15245 && alp(i)<1.13954
%         w_temp = w(i,:);
%         ind = [2 1 3 4 6 5];
%         w(i,:) = w_temp(ind);
%     end
end
w1=real(w);
w2=imag(w);
colorarry = {'b','r','k'};
colorarry2 = {'b','k','r'};

% Uncomment to plot k gap when eps_rho=0 and eps_k>0
figure
hold on
for j=N+1:2*N
    plot(alp,w1(:,j),colorarry{j-N})
%     scatter(alp,w2(:,j),20,colorarry{j-N},'filled')
end
% line([-1.31031,-1.31031],[-max(max(w2))*1.2,Omega/2],'Color','g');
% line([-1.02615,-1.02615],[-max(max(w2))*1.2,Omega/2],'Color','g');
% line([1.11435,1.11435],[-max(max(w2))*1.2,Omega/2],'Color','g','Markersize',30);
% line([1.20249,1.20249],[-max(max(w2))*1.2,Omega/2],'Color','g','Markersize',30);

% Uncomment to plot band gap when eps_rho>0 and eps_k=0
figure
hold on
alp1 = alp(1:64);
alp2 = alp(65:137);
alp3 = alp(138:end);
for j=N+1:2*N
    plot(alp1,w1(1:64,j),colorarry{j-N})
    plot(alp2,w1(65:137,j),colorarry2{j-N})
    plot(alp3,w1(138:end,j),colorarry{j-N})
%     scatter(alp,w2(:,j),20,colorarry{j-N},'filled')
end
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

% Choose a legend and adjust the x,y lim
% legend()
lgd =legend('Re(\omega)','','','','','','','','','','','','Band gap')
% % lgd =legend('Re(\omega)','Im(\omega)')
lgd.FontSize = 10;
legend('Location','southeast');
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
    [w, cnd] = hill_exp(T,M,N); % what is cnd?
    w_real=real(w);
    % w_real is a vector of size 2N=6
    % Sort w_real in order to yield the correct bands
    [w_real,ind]=sort(w_real);
    w = w(ind);
end   


