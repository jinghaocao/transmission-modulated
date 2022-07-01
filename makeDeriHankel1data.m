function deriHdata = makeDeriHankel1data(N_multi,z)
deriHdata=zeros(1,2*N_multi+1);

for n=-(N_multi):N_multi     
   deriHdata(N_multi+1+n) = 1/2*(besselh(n-1,1, z) - besselh(n+1,1, z));
end

end