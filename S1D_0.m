function sum=S1D_0(k,alp,L,v_zeta)



sum=0;


C=-psi(1);
%C=1;

b0=alp(1);
%t0=asin(b0/k);
g0=sqrt(b0^2-k^2);
    
    if  (b0^2-k^2)<0
    g0=-g0;
    
    end

    
sum=sum-1-2*1i/pi*( C +log(k*L/(4*pi))  )-2*1i/g0/L-2*1i*(k^2+2*b0^2)/L*(L/(2*pi))^3*v_zeta;


sum_temp=0;
N_lattice = 70;
for m=1:N_lattice
    
    
    bm=alp(1)+m*2*pi/L;
    bm2=alp(1)+(-m)*2*pi/L;
    
    %tm=asin(bm/k);
    %tm2=asin(bm2/k);

    
    gm=sqrt(bm^2-k^2);
    gm2=sqrt(bm2^2-k^2);
    if  (bm^2-k^2)<0
    gm=-gm;
    
    end
    if  (bm2^2-k^2)<0
    gm2=-gm2;
    
    end
    
    sum_temp=sum_temp+1/gm+1/gm2 - 2*(L/(2*pi))/m - (k^2+2*b0^2)/m^3*(L/(2*pi))^3;
    %sum_temp=sum_temp+exp(-2*1i*n*tm)/gm/L+exp(2*1i*n*tm2)/gm2/L-(-1)^n/(m*pi)*(k*L/(4*m*pi))^(2*n);
    
    
    
    
    

end

sum=sum-2*1i/L*sum_temp;




