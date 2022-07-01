function out = Mfunc(t,delta,kappab,rhob,vol,C,rhot,sqrtkappat,w3)
    Rho = diag(rhot(t));
    Rinv = diag(1./rhot(t));
    K = diag(sqrtkappat(t));
    W3 = diag(w3(t));
    out = delta*kappab/rhob*inv(diag(vol))*K*Rho*C*K*Rinv + W3;
end