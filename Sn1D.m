function sum=Sn1D(n,k,alp,L,Num_m)



sum=0;


b0=alp(1);
t0=asin(b0/k);
g0=sqrt(b0^2-k^2);
    
    if  (b0^2-k^2)<0
    g0=-g0;
    
    end

sum=sum+(-1)*2*1i*exp(-2*1i*t0)/g0/L;

sum_temp=0;
for m=1:Num_m
    
    
    bm=alp(1)+m*2*pi/L;
    bm2=alp(1)+(-m)*2*pi/L;
    
    tm=asin(bm/k);
    tm2=asin(bm2/k);

    
    gm=sqrt(bm^2-k^2);
    gm2=sqrt(bm2^2-k^2);
    if  (bm^2-k^2)<0
    gm=-gm;
    gm2=-gm2;
    end
    
    
    sum_temp=sum_temp+exp(-2*1i*n*tm)/gm/L+exp(2*1i*n*tm2)/gm2/L-(-1)^n/(m*pi)*(k*L/(4*m*pi))^(2*n);
    
    
    
    
    

end

sum=sum+(-2*1i)*sum_temp;

sum=sum+(-1)*(2*1i)*(-1)^n/pi*(L*k/(4*pi))^(2*n)*zeta(2*n+1)+1i/(n*pi);

sum_temp2=0;
for m=1:n
   
    sum_temp2=sum_temp2+(-1)^m*2^(2*m)*factorial(n+m-1)/(factorial(2*m)*factorial(n-m))*(2*pi/k/L)^(2*m)*bernoulli(2*m, alp(1)*L/(2*pi));
    
    
    
end


sum=sum+1i/pi*sum_temp2;


