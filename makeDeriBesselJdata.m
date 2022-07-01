function deriJdata = makeDeriBesselJdata(N_multi,z)
deriJdata=zeros(1,2*N_multi+1);

for n=-(N_multi):N_multi
    deriJdata(N_multi+1+n) = 1/2*(besselj(n-1, z) - besselj(n+1, z));
end
end