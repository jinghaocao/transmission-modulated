function Jdata = makeBesselJdata(N_multi,z)

Jdata=zeros(1,2*N_multi+1);

for n=-N_multi:N_multi
   
    Jdata(N_multi+1+n) = besselj(n,z);
    
end

end