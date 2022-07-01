function JHc1c2data = makeJHc1c2expdata(k0,c1,c2,N_multi)
c1c2 = c1-c2;
c2c1 = c2-c1;

normc1c2 = norm(c1c2);
thc1c2 = angle( c1c2(1) + 1i*c1c2(2) );
thc2c1 = angle( c2c1(1) + 1i*c2c1(2) );


Jc1c2data=makeBesselJdata(2*N_multi, k0*normc1c2  );
Hc1c2data=makeHankel1data(2*N_multi, k0*normc1c2  );

Jc2c1data=Jc1c2data;
Hc2c1data=Hc1c2data;

for n=-2*N_multi:2*N_multi

   Jc1c2data( n + 2*N_multi+1 )=Jc1c2data( n + 2*N_multi+1 )*exp(1i*n*thc1c2);
   Hc1c2data( n + 2*N_multi+1 )=Hc1c2data( n + 2*N_multi+1 )*exp(1i*n*thc1c2);
   
   Jc2c1data( n + 2*N_multi+1 )=Jc2c1data( n + 2*N_multi+1 )*exp(1i*n*thc2c1);
   Hc2c1data( n + 2*N_multi+1 )=Hc2c1data( n + 2*N_multi+1 )*exp(1i*n*thc2c1);
   
end

JHc1c2data = [Jc1c2data;Hc1c2data;Jc2c1data;Hc2c1data];
