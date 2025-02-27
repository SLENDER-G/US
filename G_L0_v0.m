function [x, objective,criterion] = G_L0_v0(Phi,DAS_RF, b, eps, mu, lambda, beta1, beta2, beta3, maxiter)
global b u_k lambd_k1 lambd_k2 lambd_k3 
x_k = zeros(484137,1);
u_k = zeros(484137,1);
X = zeros(1251,387);
lambd_k1 = zeros(1251,387);
lambd_k2 = zeros(484137,1);
lambd_k3 = zeros(484137,1);
objective = 1e6;
options = [];
options.display = 'none';
options.maxFunEvals = 100;
options.Method = 'lbfgs';
for i=1:maxiter
    %%%%%%%%%%%%%%%%%%%%%%%%%% update u %%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_k = (DAS_RF+beta1*(X(:)+lambd_k1(:)/beta1)+beta2*(u_k-lambd_k2/beta2))/(1+beta1+beta2);
   % x_k = minFunc(@(x) L0_FunCtion(b, Phi, x, u_k, lambd_k1(:), lambd_k2, beta1, beta2, X(:)), zeros(484137,1), options);
    %%%%%%%%%%%%%%%%%%%%%%%%%% update G %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % argmin_X lambda||X||+ beta1/2*||X - x + lambd_k1/beta1||^2
    Xk = reshape(x_k,[1251 387])-lambd_k1/beta1;
    X = nuclear_norm_matrix(Xk,lambda/beta1,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%% update v %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %cof_A = beta3*abs(u_k).*abs(u_k);
    cof_A = beta3*u_k.*u_k;
    cof_B = mu-lambd_k3.*abs(u_k);
    v_k=boxproj(cof_B./cof_A);
    %%%%%%%%%%%%%%%%%%%%%%%%%% update u_k %%%%%%%%%%%%%%%%%%%%%%%%%%%
    q = x_k+lambd_k2/beta2;
    p = (2.*abs(q)-v_k.*lambd_k3)./(2+beta3*v_k.*v_k);
    u_k=threadholding_l1_w(q,p);
    %%%%%%%%%%%%%%%%%%%%%%%%%% update lambd %%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambd_k1 = lambd_k1+beta1*(X-reshape(x_k,[1251 387]));
    lambd_k2 = lambd_k2+beta2*(x_k-u_k);
    lambd_k3 = lambd_k3+beta3*(v_k.*abs(u_k));
    %%%%%%%%%%%%%%%%%%%% stopping criterion %%%%%%%%%%%%%%%%%%%%%%%
    objective(i+1) = 0.5*norm(b-full(Phi*sparse(double(x_k))))^2 +mu * sum(1 - v_k)+lambda*sum(abs(X(:)))...
        +(beta1/2)*norm(X(:)-x_k+(lambd_k1(:)/beta1))^2 ...
        +(beta2/2)*norm(x_k-u_k+(lambd_k2/beta2))^2 ...
        +(beta3/2)*norm(v_k.*u_k+(lambd_k3/beta3))^2;
    criterion(i) = abs(objective(i+1)-objective(i))/objective(i);
    disp(['iteration ',num2str(i),', criterion: ',num2str(criterion(i))])
    
    if ( criterion(i) < eps )
        x(:,i) = x_k;
        break
    end
        x(:,i) = x_k;
end
objective = objective(2:end);

end

