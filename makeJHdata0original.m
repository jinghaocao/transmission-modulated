function JHdata = makeJHdata0original(k0,R,N_multi)

%%% collect values of the bessel and hankel functions
Jdata_k0R=makeBesselJdata(N_multi, k0*R  );
Hdata_k0R=makeHankel1data(N_multi, k0*R  );

JHdata = zeros(2,2*N_multi+1);
JHdata(1,:) = Jdata_k0R;
JHdata(2,:) = Hdata_k0R;