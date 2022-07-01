function JHdata = makeJHdata(k0,kb,R,N_multi)


%%% collect values of the bessel and hankel functions
Jdata_k0R=makeBesselJdata(N_multi, k0*R  );
Jdata_kbR=makeBesselJdata(N_multi, kb*R  );
Hdata_k0R=makeHankel1data(N_multi, k0*R  );
Hdata_kbR=makeHankel1data(N_multi, kb*R  );
dJdata_k0R=makeDeriBesselJdata(N_multi,k0*R);
dJdata_kbR=makeDeriBesselJdata(N_multi,kb*R);
dHdata_k0R=makeDeriHankel1data(N_multi,k0*R);

JHdata = zeros(7,2*N_multi+1);
JHdata(1,:) = Jdata_k0R;
JHdata(2,:) = Jdata_kbR;
JHdata(3,:) = Hdata_k0R;
JHdata(4,:) = Hdata_kbR;
JHdata(5,:) = dJdata_k0R;
JHdata(6,:) = dJdata_kbR;
JHdata(7,:) = dHdata_k0R;