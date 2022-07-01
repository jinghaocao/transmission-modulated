function Sn=lattice_Sn(n,k,alp,L1x,L2,N_lattice,d_zeta)


    
    
if mod(n,2) == 0
    %number is even
      
    if n==0
      
      Sn=S1D_0(k,alp,L1x,d_zeta(1))+DeltaSn(0,k,alp,L1x,L2,N_lattice);  
    else
      l=n/2;
      Sn=S1D_2n(l,k,alp,L1x,d_zeta(l))+DeltaSn(n,k,alp,L1x,L2,N_lattice);
      %Sn=Sn*(-1)^n;
    end
else
  %number is odd
  l=(n+1)/2;
  Sn=S1D_2n_1(l,k,alp,L1x,d_zeta(l))+DeltaSn(n,k,alp,L1x,L2,N_lattice);
  %Sn=Sn*(-1)^n;
end




