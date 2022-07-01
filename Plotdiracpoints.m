D =1;
%%% Set reciprocal vectors
x = sqrt(3);
b1 = 4*pi/(3*D)*[x/2, -1/2];
b2 = 4*pi/(3*D)*[0, 1];

%%% Set symmetry points in reciprocal space
pointG = [0,0];
pointM = 1/2*b1;
pointK = 2/3*b1 + 1/3*b2;

Omega = 0.2;
epsr = 0.5;
steps = 50;
alps_dirac =[];
w_real_dirac = [];
R = @(t) [cos(t) -sin(t); sin(t) cos(t)];
% points0 = [pointK,0.5*(pointK+pointM),0.5*(pointG+pointK)];
for t = 0:5
    pointKi= (R(2*pi/6*t)*pointK')';
    pointMi = (R(2*pi/6*t)*0.5*(pointK+pointM)')';
    pointGi = (R(2*pi/6*t)*0.5*(pointK+pointG)')';
    [lengthAB1,alps_real,w_real,w_imag] = HonelyComb_Banddata(pointMi,pointKi,Omega,epsr,steps);
    w_real_dirac = [w_real_dirac, w_real];
    alps_dirac = [alps_dirac,(0:steps-1)*lengthAB1/steps+3*t];
    [lengthAB2,alps_real,w_real,w_imag] = HonelyComb_Banddata(pointKi,pointGi,Omega,epsr,steps);
    alps_dirac = [alps_dirac,(0:steps-1)*lengthAB2/steps+lengthAB1+3*t];
    w_real_dirac = [w_real_dirac, w_real];
end

fsz = 16;
lfsz = 14;
alw = 0.75;
lw = 1.2;

figure 
N = 6;
hold on
for I_w = 1:N*2
    plot(alps_dirac,w_real_dirac(I_w,:),'.b','LineWidth',lw)
end