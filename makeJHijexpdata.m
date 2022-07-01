function JHijdata = makeJHijexpdata(k0,c,N_multi)
N = size(c,1);
NN = N*(N-1)/2; % Number of pairs of particles
NN_multi = 4*N_multi + 1;
JHijdata = zeros(2,NN_multi,N,N);
for i = 1:N
    for j = i+1:N
        cicj = c(i,:)-c(j,:);

        normcicj = norm(cicj);
        thcicj = angle( cicj(1) + 1i*cicj(2) );

        Jcicjdata=makeBesselJdata(2*N_multi, k0*normcicj  );
        Hcicjdata=makeHankel1data(2*N_multi, k0*normcicj  );
        
        for n=-2*N_multi:2*N_multi

           Jcicjdata( n + 2*N_multi+1 )=Jcicjdata( n + 2*N_multi+1 )*exp(1i*n*thcicj);
           Hcicjdata( n + 2*N_multi+1 )=Hcicjdata( n + 2*N_multi+1 )*exp(1i*n*thcicj);

        end

        JHijdata(:,:,i,j) = [Jcicjdata;Hcicjdata];
    end
end

