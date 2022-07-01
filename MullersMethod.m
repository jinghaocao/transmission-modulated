function [root, fRoot, iIter, iterRoots] = MullersMethod(f, z0, z1, z2, iterMax, distTol, fTol)
% MullersMethod
%
% Overview:
%   Calculates the root of a complex valued function.
%
% Input:
%   f:          A function for which we want to determine a complex root.
%   z0:         A complex number. initial guess at the root.
%   z1:         An initial guess at the root.
%   z2:         An initial guess at the root.
%   iterMax:    The maxmimum number of iterations to perform.
%   distTol:    When the distance between the values of the roots in
%               successive successive iterations falls below this value the
%               process is terminated.
%   fTol:       When the value of the function evaluated at an approximate
%               root falls below this value the process is terminated.
%
% Output:
%   root:       The approximate complex root of f.
%   fRoot:      The value of the function at the root.
%   iterN:      The numer of iterations used to obtain the root.
%   iterRoots:  A vector containing the approximate roots obtained on each
%               iteration.
%
% References:
%   Mathematical and Computational Methods in Photonics - Tutorial Notes
%
% Authors:
%   Habib Ammari, Brian Fitzpatrick, Matias Ruiz, Sanghyeon Yu.

% --------------------------------------------------
% Ensure the initial guesses are distinct.
% --------------------------------------------------
if z0==z1 || z1==z2 || z2==z0
    error('The initial points must be distinct.')
end
    
% --------------------------------------------------
% Define the divided differences for f.
% --------------------------------------------------
f0 = feval(f, z0);
f1 = feval(f, z1);
f2 = feval(f, z2);

iterRoots = zeros(iterMax, 1);

% --------------------------------------------------
% Iterate until we find a root or reach the maxmimum number
% of iterations.
% --------------------------------------------------
for iIter = 1:iterMax
    q = (z2 - z1)/(z1 - z0);
    a = q*f2 - q*(1+q)*f1 + q^2*f0;
    b = (2*q + 1)*f2 - (1 + q)^2*f1 + q^2*f0;
    c = (1 + q)*f2;
    
    b24mac = b^2 - 4*a*c;

    denom1 = b + sqrt(b24mac);
    denom2 = b - sqrt(b24mac);
    
    if abs(denom1) > abs(denom2)
        z3 = z2 - (z2 - z1)*(2*c/denom1);
    else
        z3 = z2 - (z2 - z1)*(2*c/denom2);
    end
    
    iterRoots(iIter) = z3;

    f3 = feval(f, z3);
 
    if (abs(z3 - z2) < distTol && abs(f3) < fTol )
        root = z3;
        fRoot = f3;
        iterRoots(iIter:end) = [];      
        return
    end     
    
    z0 = z1;
    z1 = z2;
    z2 = z3;
    
    f0 = f1;
    f1 = f2;
    f2 = f3;    
end
root = z3;
warning('No root was found before the maximumn number of iterations was reached.')