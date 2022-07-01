%%% Make multipole approximation of single-layer (part of operator A) 

function matS = makeS_1D(k0,R,alp,L,d_zeta,Jdata,Hdata,JHijdata,N,N_multi)
NN = 2*N_multi + 1;
matS = zeros(NN*N);
%%% Collect lattice sum data
data_Sn_posi = zeros(2*N_multi+1,1);
data_Sn_nega = zeros(2*N_multi,1);
for j=1:2*N_multi+1    
data_Sn_posi(j)=lattice_Sn_1D(j-1,k0,alp,L,d_zeta);
end
for j=1:2*N_multi   
data_Sn_nega(j)=lattice_Sn_1D(j,k0,-alp,L,d_zeta);
end

for i=1:N
    %%% --------Create single-particle operator (matrix Sk0)--------
    Sk0=zeros(2*N_multi+1);
    %%% define constant c
    const_i=(-1/2)*1i*pi*R(i);
    for n=-N_multi:N_multi
       Jk0= Jdata(i,N_multi+1+n);
       Hk0= Hdata(i,N_multi+1+n);
       Sk0(n+N_multi+1,n+N_multi+1)=const_i*Jk0*Hk0;
       for m=-N_multi:N_multi
            if (n-m) >= 0        
                Snm = data_Sn_posi(n-m+1);
            else
                jj=-(n-m);
                Snm = data_Sn_nega(jj);
            end             
            Jk0m = Jdata(i,N_multi+1+m);        
            Sk0(m+N_multi+1,n+N_multi+1)=Sk0(m+N_multi+1,n+N_multi+1) + const_i*Jk0*(-1)^(n-m)*Snm*Jk0m;               
       end   
    end
    
%%% --------Create particle-to-particle operators--------
    sgn = (-1).^(-2*N_multi:2*N_multi);

    i_I = NN*(i-1)+1:NN*i;
    for j = 1:N
        const_j = (-1/2)*1i*pi*R(j);
        if i == j      
            matS(i_I,i_I) = Sk0;
        else
            j_I = NN*(j-1)+1:NN*j;
            S_ij=zeros(NN);
            if j > i
                Jcicjdata = JHijdata(1,:,i,j);
                Hcicjdata = JHijdata(2,:,i,j);
            else 
                Jcicjdata = JHijdata(1,:,j,i).*sgn;
                Hcicjdata = JHijdata(2,:,j,i).*sgn;
            end

            for n=-N_multi:N_multi
               Jk0= Jdata(j,N_multi+1+n);
               for m=-N_multi:N_multi
                   Jk0m = Jdata(i,N_multi+1+m);

                   sumQcicj = 0;
                   for l = -N_multi:N_multi
                       if (n-l) >= 0
                            Snl = data_Sn_posi(n-l+1);
                       else
                            jj=-(n-l);
                            Snl = data_Sn_nega(jj);
                       end
                       sumQcicj = sumQcicj + (-1)^(n-l)*Snl*(-1)^(l-m)*Jcicjdata( l-m + 2*N_multi+1 );
                   end

                   tempDiDj = Hcicjdata(n-m + 2*N_multi+1) + sumQcicj;

                   S_ij(m+N_multi+1,n+N_multi+1)=const_j*Jk0*Jk0m*tempDiDj;       
               end
            end
            matS(i_I,j_I) = S_ij; % Make sure transpose is correct
        end
    end
end