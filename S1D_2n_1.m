function sum=S1D_2n_1(n,k,alp,L,v_zeta)


%load b_data
sum=0;


b0=alp(1);
t0=asin(b0/k);
g0=sqrt(b0^2-k^2);
    
    if  (b0^2-k^2)<0
    g0=-g0;
    
    end

sum=sum+2*1i*exp(-1i*(2*n-1)*t0)/g0/L;

sum_temp=0;
N_lattice = 70;
for m=1:N_lattice
    
    
    bm=alp(1)+m*2*pi/L;
    bm2=alp(1)+(-m)*2*pi/L;
    
    tm=asin(bm/k);
    tm2=asin(bm2/k);

    
    gm=sqrt(bm^2-k^2);
    gm2=sqrt(bm2^2-k^2);
    if  (bm^2-k^2)<0
    gm=-gm;
    
    end
    if  (bm2^2-k^2)<0
    gm2=-gm2;
    
    end
    
    sum_temp=sum_temp+exp(-1i*(2*n-1)*tm)/gm/L-exp(1i*(2*n-1)*tm2)/gm2/L+1i*(-1)^n*b0*L*n/(m*pi)^2*(k*L/(4*m*pi))^(2*n-1);
    
    
    
    
    

end

sum=sum+2*1i*sum_temp;

sum=sum+2*(-1)^n*b0*L*n/pi^2*(L*k/(4*pi))^(2*n-1)*v_zeta;


sum_temp2=0;


for m=0:n-1
    %%%sum_temp2=sum_temp2+(-1)^m*2^(2*m)*factorial(n+m-1)/(factorial(2*m+1)*factorial(n-m-1))*(2*pi/k/L)^(2*m+1)*bernoulli(2*m+1, alp(1)*L/(2*pi));
    %%%sum_temp2=sum_temp2+(-1)^m*2^(2*m)*gamma(n+m)/(gamma(2*m+2)*gamma(n-m))*(2*pi/k/L)^(2*m+1)*b_data(2*m+1);
    %%%sum_temp2=sum_temp2+(-1)^m*2^(2*m)*factorial(n+m-1)/(factorial(2*m+1)*factorial(n-m-1))*(2*pi/k/L)^(2*m+1)*bernoulli_poly2(2*m+1, alp(1)*L/(2*pi));
    
    
    sum_temp2=sum_temp2+(-1)^m*2^(2*m)*exp(gammaln(n+m-1 +1))/( exp(gammaln(2*m+1 +1))*exp(gammaln(n-m-1 +1)))*(2*pi/k/L)^(2*m+1)*bernoulli_poly2(2*m+1, alp(1)*L/(2*pi));
    %sum_temp2=sum_temp2+(-1)^m*2^(2*m)*exp(gammaln(n+m-1 +1))/( exp(gammaln(2*m+1 +1))*exp(gammaln(n-m-1 +1)))*(2*pi/k/L)^(2*m+1)*1;
    
    
    %%%temp1=bernoulli(2*m+1, alp(1)*L/(2*pi))
    %%%temp2=b_data(2*m+1)
    
end


sum=sum-2/pi*sum_temp2;


