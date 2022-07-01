%%% Make multipole approximation of single-layer (part of operator A) 

function matS = makeS(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multi,N_lattice)
NN = 2*N_multi + 1;
matS = zeros(NN*N);

%%% Collect lattice sum data
data_Sn_posi = zeros(2*N_multi+1,1);
data_Sn_nega = zeros(2*N_multi,1);
for j=1:2*N_multi+1    
    data_Sn_posi(j)=lattice_Sn(j-1,k0,alp,L1x,L2,N_lattice,d_zeta);
end
for j=1:2*N_multi   
    data_Sn_nega(j)=lattice_Sn(j,k0,[-alp(1), alp(2)],L1x,L2,N_lattice,d_zeta);
end

%%% define constant c
const=(-1/2)*1i*pi*R;

%%% --------Create single-particle operator (matrix Sk0)--------
Sk0=zeros(2*N_multi+1);

%%% Collect bessel and hankel functions data
Jdata_k0R = JHdata(1,:);
Hdata_k0R = JHdata(2,:);

for n=-N_multi:N_multi
   Jk0= Jdata_k0R(N_multi+1+n);
   Hk0= Hdata_k0R(N_multi+1+n);
       
   Sk0(n+N_multi+1,n+N_multi+1)=const*Jk0*Hk0;
       
   for m=-N_multi:N_multi
        if (n-m) >= 0        
            Snm = data_Sn_posi(n-m+1);
        else
            jj=-(n-m);
            Snm = data_Sn_nega(jj);
        end             
        Jk0m = Jdata_k0R(N_multi+1+m);        
        Sk0(m+N_multi+1,n+N_multi+1)=Sk0(m+N_multi+1,n+N_multi+1) + const*Jk0*(-1)^(n-m)*Snm*Jk0m;               
   end   
end

%%% --------Create particle-to-particle operators--------
sgn = (-1).^(-2*N_multi:2*N_multi);
for i = 1:N
    i_I = NN*(i-1)+1:NN*i;
    matS(i_I,i_I) = Sk0;
    for j = i+1:N
        j_I = NN*(j-1)+1:NN*j;
        
        S_DiDj=zeros(NN);
        S_DjDi=zeros(NN);
               
        Jcicjdata = JHijdata(1,:,i,j);
        Hcicjdata = JHijdata(2,:,i,j);
        Jcjcidata = Jcicjdata.*sgn;
        Hcjcidata = Hcicjdata.*sgn;
        
        for n=-N_multi:N_multi
           Jk0= Jdata_k0R(N_multi+1+n);
           for m=-N_multi:N_multi
               Jk0m = Jdata_k0R(N_multi+1+m);

               sumQcicj = 0;
               sumQcjci = 0;
               for l = -N_multi:N_multi
                   if (n-l) >= 0
                        Snl = data_Sn_posi(n-l+1);
                   else
                        jj=-(n-l);
                        Snl = data_Sn_nega(jj);
                   end
                   sumQcicj = sumQcicj + (-1)^(n-l)*Snl*(-1)^(l-m)*Jcicjdata( l-m + 2*N_multi+1 );
                   sumQcjci = sumQcjci + (-1)^(n-l)*Snl*(-1)^(l-m)*Jcjcidata( l-m + 2*N_multi+1 );
               end

               tempDiDj = (-1)^(n-m)*Hcicjdata(n-m + 2*N_multi+1) + sumQcicj;
               tempDjDi = (-1)^(n-m)*Hcjcidata(n-m + 2*N_multi+1) + sumQcjci;

               S_DiDj(m+N_multi+1,n+N_multi+1)=const*Jk0*Jk0m*tempDiDj;       
               S_DjDi(m+N_multi+1,n+N_multi+1)=const*Jk0*Jk0m*tempDjDi;
           end
        end
        matS(i_I,j_I) = S_DjDi;
        matS(j_I,i_I) = S_DiDj;
    end
end