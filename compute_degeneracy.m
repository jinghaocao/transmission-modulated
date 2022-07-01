%% calculate the degeneracy between the the first and second band

function alpha_wanted = compute_degeneracy(c,R,L,Omega,kappa_b,rho_b,vol,delta) 
    T = 2*pi/Omega;    
    N = size(R,1);
    eigenvalues = zeros(6,1);
    N_multipole=4;
    error = 0.01;
    d_zeta=makezetadata;
    k0 = 0.00001;
    [Jdata,Hdata] = makeJHdata0(k0,R,N_multipole);
    JHijdata = makeJHijexpdata(k0,c,N_multipole);
    alpha_wanted = 0;
    for alpha = -2.3:0.0005:-2.2
        C = makeC_1D(k0,R,alpha,L,d_zeta,Jdata,Hdata,JHijdata,N,N_multipole);
        M = @(t) delta*kappa_b/rho_b*inv(diag(vol))*C;
        [eigenvalues, cnd] = hill_exp(T,M,N);
        eigenvalues = sort(real(eigenvalues));
        if abs(eigenvalues(1)-eigenvalues(2))< error
            alpha_wanted = alpha;
            eigenvalues_wanted = eigenvalues;
            error = abs(eigenvalues_wanted(1)-eigenvalues_wanted(2));
        end
    end
end
