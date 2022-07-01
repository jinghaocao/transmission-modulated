function [Jdata,Hdata] = makeJHdata0(k0,R,N_multi)

%%% collect values of the bessel and hankel functions
N = size(R,1);
Jdata = zeros(N,2*N_multi+1);
Hdata = zeros(N,2*N_multi+1);
for i = 1:N
    Jdata(i,:)= makeBesselJdata(N_multi, k0*R(i));
    Hdata(i,:)= makeHankel1data(N_multi, k0*R(i));
end