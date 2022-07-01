%%% Make capacitance matrix

function matC = makeC_1D(k0,R,alp,L,d_zeta,Jdata,Hdata,JHijdata,N,N_multi)
NN = 2*N_multi + 1;
M = NN*N;
matS = makeS_1D(k0,R,alp,L,d_zeta,Jdata,Hdata,JHijdata,N,N_multi);
matC = zeros(N,N);
for j = 1:N
    phi_j = zeros(M,1); 
    phi_j(NN*(j-1)+N_multi+1) = 1;
    psi_j = matS\phi_j;
    for i = j:N
        phi_i = zeros(M,1); 
        phi_i(NN*(i-1)+N_multi+1) = 1;
        matC(i,j) = -2*pi*R(i)*phi_i'*psi_j;
        matC(j,i) = conj(matC(i,j));
    end
end