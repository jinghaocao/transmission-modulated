function W = bloch_mode(M,M0,N,T)
    W = zeros(2*N,2*N);
    Z = zeros(N,N);
    I = eye(N,N);
    MM = @(t,y) [Z, I; -M(t), Z]*y;
    [Eig_0,D_0] = eig([Z, I; -M0(0), Z]);
    
    
    % Thea's initial condition
    II = eye(2*N,2*N);
    X = zeros(2*N,2*N);
    for j = 1 : 2*N
       [~,y] = ode45(MM,[0,T],II(:,j)); 
       X(:,j) = y(end,:)'; 
    end
    [V,d] = eig(X,'vector');
    f_exp = log(d)/1i/T;
    [f_exp_real, ind] = sort(real(f_exp),'descend');
    f_exp = f_exp(ind);
    V = V(:, ind);
    d = d(ind);
    Eig_0 = V;
    TT = 5*T;
    figure
    hold on
    for j = 1:N*2
        y_0 = Eig_0(:,j);
        [~,yI] = ode45(MM,[0,TT],y_0);
        W(:,j)=yI(end,:).';
%         NN = length(yI);
%         y = zeros(1,NN);
%         for i  = 1:NN
%             y(i) = norm(yI(i,1:6));
%         end 
%         plot([1:NN],y);
    end
end