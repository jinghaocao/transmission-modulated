function sum=DeltaSn(n,omega,alp,L1x,L2,Num_m)


L2x=L2(1); 
L2y=L2(2); 

sum=0;




for m=-Num_m:Num_m
    
    
    bm=alp(1)+m*2*pi/L1x;
 
    
    gm=sqrt(omega^2-bm^2);
    
    expitm = (gm+1i*bm)/omega;
    
    mu = exp(1i*(alp*L2'));
    Pm = exp(1i*gm*L2y);
    Qm = exp(1i*bm*L2x);
    
    
    sum=sum+1/gm*(  expitm^n/(1/mu*Qm/Pm - 1)  + (-1)^n*expitm^(-n)/( mu/Qm/Pm - 1)  );
    
    
 
end

sum=sum*2/L1x;


