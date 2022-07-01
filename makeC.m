%%% Make capacitance matrix

function matC = makeC(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multi,N_lattice)
NN = 2*N_multi + 1;
M = NN*N;
matS = makeS(k0,R,alp,L1x,L2,d_zeta,JHdata,JHijdata,N,N_multi,N_lattice);
matC = zeros(N,N);
for j = 1:N
    phi_j = zeros(M,1); 
    phi_j(NN*(j-1)+N_multi+1) = 1;
    psi_j = matS\phi_j;
    for i = j:N
        phi_i = zeros(M,1); 
        phi_i(NN*(i-1)+N_multi+1) = 1;
        matC(i,j) = -2*pi*R*phi_i'*psi_j;
        matC(j,i) = conj(matC(i,j));
    end
end
