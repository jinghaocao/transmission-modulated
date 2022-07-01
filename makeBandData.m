%%% Plot the band gap
function [w_real,w_imag,alps]= makeBandData(pointM,c,Omega,delta,kappa_b,rho_b,R,L,vol,rhot,sqrtkappat,w3,N_a)
    %%% Set the discretization or truncation parameters
    T = 2*pi/Omega;
    N = size(R,1);
    N_multipole=4;
    % alps = linspace(-pointM,pointM,N_a);
    alps = linspace(-pointM,pointM,N_a);

    %%% Precompute values of zeta function (used in computation of lattice sums, just for acceleration)
    d_zeta=makezetadata;
    % BandData = zeros(2*N,N_a);
    % cnd = zeros(1,N_a);

    k0 = 0.00001;
    [Jdata,Hdata] = makeJHdata0(k0,R,N_multipole);
    JHijdata = makeJHijexpdata(k0,c,N_multipole);
    
    parfor alp_ind = 1:N_a
        %%% Set quasi-periodic parameter \alpha
        alp =  alps(alp_ind);
        %%% Compute the operator A
        C = makeC_1D(k0,R,alp,L,d_zeta,Jdata,Hdata,JHijdata,N,N_multipole);
        M = @(t) Mfunc(t,delta,kappa_b,rho_b,vol,C,rhot,sqrtkappat,w3);
        [BandData(:, alp_ind), cnd(alp_ind)] = hill_exp(T,M,N);
    end  

    %% Plots the solution and saves good images
    fsz = 16;
    lfsz = 14;
    alw = 0.75;
    lw = 1.2;

    [w_real_temp, I_sort] = sort(real(BandData)); % sort(real(bands),1);
    for I_a=1:N_a
        for I_w=1:2*N
            w_imag_temp(I_w,I_a) = imag(BandData(I_sort(I_w,I_a),I_a)); % sort(imag(bands),1);
        end
    end

    Ias = find(max(abs(w_imag_temp))>1e-4);
    [w_i_sort, I_sort] = sort(w_imag_temp(N+1:2*N,Ias));

    w_imag = sort(w_imag_temp);
    w_real = w_real_temp;

end