function Hdata = makeHankel1data(N_multi,z)
Hdata=zeros(1,2*N_multi+1);

for n=-N_multi:N_multi
    Hdata(N_multi+1+n) = besselh(n,1,z);
end

end