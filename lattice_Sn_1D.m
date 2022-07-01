function Sn=lattice_Sn_1D(n,k,alp,L1x,d_zeta)


    
    
if mod(n,2) == 0
    %number is even
      
    if n==0
      
      Sn=S1D_0(k,alp,L1x,d_zeta(1));
    else
      l=n/2;
      Sn=S1D_2n(l,k,alp,L1x,d_zeta(l));
      %Sn=Sn*(-1)^n;
    end
else
  %number is odd
  l=(n+1)/2;
  Sn=S1D_2n_1(l,k,alp,L1x,d_zeta(l));
  %Sn=Sn*(-1)^n;
end




