%%% Copyright (C) 2001 - 2002 - R.Coisson & G. Vernizzi, Parma University
% Copyright (C) X. K. Yang
% Copyright (C) 2010 - 2012 - N. O. Strelkov, NRU MPEI
% Copyright (C) 2011 - DIGITEO - Michael Baudin
% This file must be used under the terms of the GPL Version 2
% http://www.gnu.org/licenses/gpl-2.0.html
	
function [nu1, nu2] = mathieuexp(a,q,ndet)
% [nu,c]=mathieuexp(a,q,matrix_dimension)
% This program evaluates the characteristic exponent nu, 
% corresponding to solutions of Mathieu Equation 
% y''(t)+(a-2q cos(2t)) y(t)=0;
% where a and q  are fixed real (?) variables.
%
% ndet is a positive integer number: it is the matrix dimension used 
% in the algorithm. Precision increases by increasing ndet. 
% Typical value is ndet=12
% 
% The alghoritm consider two different cases: 
% a=(2k)^2 or not (k integer).
% nu is such that its real part belongs to the interval [0,2[ ([0,1[?)  
% Of course, every other solutions are obtained by the formula
% +-nu+ 2*k, with k integer.
% 
% by R.Coisson and G.Vernizzi, Parma University, 2000-2001.
% 
% 
% if (sqrt(abs(a))/2-fix(sqrt(abs(a))/2))==0 &a>=0, %two cases considered 
% x=1;%a=(2k)^2
% else 
% x=0;%a~=(2k)^2
% end
	
    x=(a>=0)& (mod(sqrt(abs(a))/2,1)==0);
    N=2*ndet+1; % matrix dimension
    d=q./((2.*(-ndet:ndet)+x).^2-a); % off-diagonal elements 
    m=eye(N,N)+diag(d(1:N-1),1)+diag(d(2:N),-1); % matrix definition
    delta = det(m); % evaluation of determinant
    if x==1 % calculation of alpha in two cases  
        alpha = acos(2*delta-1)/pi;
    else
        alpha = 2*asin(sqrt(delta)*(sin(sqrt(a)*pi/2)))/pi;
    end
    nu1=mod(real(alpha),2)+1i*imag(alpha);
    nu2=mod(real(-alpha),2)+1i*imag(-alpha);
end


%     %PRIMA OPZIONE
%     %nu=abs(mod(real(alpha),2))+1i*imag(alpha); %just modular reduction to the solution [0,2[
%     %no change of the sign
% 
%     %SECONDA OPZIONE
%     %the same as above, but with change an overall sign such that the  imaginary part is always positive. 
%     nu=alpha*(2*(imag(alpha)>=0)-1);
%     nu=mod(real(nu),2)+1i*imag(nu);
% 
%     %TERZA OPZIONE
%     %mra=mod(real(alpha),2);
%     %nu=(1-abs(mra-1))+1i*imag(alpha)*(2*(mra<=1)-1);    %modular reduction to the solution [0,1]
%     % (with overall change sign, if necessary)
% 
% 
%     dd = (nu+2*(-ndet:ndet)).^2; %building matrix H_nu-a
%     dd_sub = q*ones(N-1,1);
%     H_nu = diag(dd-a)+diag(dd_sub,1)+diag(dd_sub,-1);
%     [U,D] = spec(H_nu); %diagonalizing
%     [minim,pos] = min(abs(diag(D))); %select eigenvector corresp. to eigenvalue=0
%     c = U(:,pos);
% 
% 
%     % for h=1:100; q(h)=(h-1)/10; a(h)=mathieu_mathieuexp(2,q(h),12); end
