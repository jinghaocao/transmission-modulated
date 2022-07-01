%% Multipole discretrization of integral operator A

function matA = makeA(k0,kb,~,alp,delta,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multi,N_lattice)
NN = 2*N_multi + 1;
matA = zeros(2*NN*N);

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
const = 1; %const=(-1/2)*1i*pi*R;

%%% --------Create single-particle operator (matrix M)--------
Sk0=zeros(2*N_multi+1);
Skb=zeros(2*N_multi+1);
dSk0=zeros(2*N_multi+1);
dSkb=zeros(2*N_multi+1);

%%% Collect bessel and hankel functions data
Jdata_k0R = JHdata(1,:);
Jdata_kbR = JHdata(2,:);
Hdata_k0R = JHdata(3,:);
Hdata_kbR = JHdata(4,:);
dJdata_k0R = JHdata(5,:);
dJdata_kbR = JHdata(6,:);
dHdata_k0R = JHdata(7,:);

for n=-N_multi:N_multi
   Jk0= Jdata_k0R(N_multi+1+n);
   Jkb= Jdata_kbR(N_multi+1+n);
   Hk0= Hdata_k0R(N_multi+1+n);
   Hkb= Hdata_kbR(N_multi+1+n);
   dHk0= dHdata_k0R(N_multi+1+n);
   dJkb= dJdata_kbR(N_multi+1+n);
       
   Sk0(n+N_multi+1,n+N_multi+1)=const*Jk0*Hk0;
   dSk0(n+N_multi+1,n+N_multi+1)=const*k0*Jk0*dHk0;
   Skb(n+N_multi+1,n+N_multi+1)=const*Jkb*Hkb;
   dSkb(n+N_multi+1,n+N_multi+1)=const*kb*Hkb*dJkb;
       
   for m=-N_multi:N_multi
      
        if (n-m) >= 0
            Snm = data_Sn_posi(n-m+1);
        else
            jj=-(n-m);
            Snm = data_Sn_nega(jj);
        end
             
        Jk0m = Jdata_k0R(N_multi+1+m);
        dJk0m= dJdata_k0R(N_multi+1+m);
        
        Sk0(m+N_multi+1,n+N_multi+1)=Sk0(m+N_multi+1,n+N_multi+1) + const*Jk0*(-1)^(n-m)*Snm*Jk0m;
        dSk0(m+N_multi+1,n+N_multi+1)=dSk0(m+N_multi+1,n+N_multi+1) + const*k0*Jk0*(-1)^(n-m)*Snm*dJk0m;
               
   end
   
end

matM=[Sk0, -Skb; delta*dSk0, -dSkb];

%%% --------Create particle-to-particle operator (matrices Lij and Lji)--------
sgn = (-1).^(-2*N_multi:2*N_multi);

for i = 1:N
    i_I = 2*NN*(i-1)+1:2*NN*i;
    matA(i_I,i_I) = matM;
    for j = i+1:N
        j_I =  2*NN*(j-1)+1:2*NN*j;
        
        S_DiDj=zeros(2*N_multi+1);
        dS_DiDj=zeros(2*N_multi+1);
        S_DjDi=zeros(2*N_multi+1);
        dS_DjDi=zeros(2*N_multi+1);

        Jcicjdata = JHijdata(1,:,i,j);
        Hcicjdata = JHijdata(2,:,i,j);
        Jcjcidata = Jcicjdata.*sgn;
        Hcjcidata = Hcicjdata.*sgn;
        for n=-N_multi:N_multi
           Jk0= Jdata_k0R(N_multi+1+n);
           for m=-N_multi:N_multi

               Jk0m = Jdata_k0R(N_multi+1+m);
               dJk0m= dJdata_k0R(N_multi+1+m);

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

               tempDiDj = (-1)^(n-m)*Hcicjdata(n-m + 2*N_multi+1)   + sumQcicj;
               tempDjDi = (-1)^(n-m)*Hcjcidata(n-m + 2*N_multi+1)   + sumQcjci;

               S_DiDj(m+N_multi+1,n+N_multi+1)=const*Jk0*Jk0m*tempDiDj;
               dS_DiDj(m+N_multi+1,n+N_multi+1)=const*k0*Jk0*dJk0m*tempDiDj;

               S_DjDi(m+N_multi+1,n+N_multi+1)=const*Jk0*Jk0m*tempDjDi;
               dS_DjDi(m+N_multi+1,n+N_multi+1)=const*k0*Jk0*dJk0m*tempDjDi;

           end

        end
        matLij=[S_DiDj, zeros(2*N_multi+1) ; delta*dS_DiDj, zeros(2*N_multi+1)];
        matLji=[S_DjDi, zeros(2*N_multi+1) ; delta*dS_DjDi, zeros(2*N_multi+1)];

        matA(i_I,j_I) = matLji;
        matA(j_I,i_I) = matLij;
    end
end