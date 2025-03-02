function [u, v, objective] = G_SR_v0(Phi,DAS_RF, b, eps, mu, lambda, beta1, beta2, maxiter)
global b v_k lambd_k1 lambd_k2  
%% loading dataset
%-- load scan and dataset
scan = linear_scan();
scan.read_file(path_scan);
dataset = us_dataset();
dataset.read_file(path_dataset);
%-- define scan based on time axis
time = (150:1400).'/dataset.sampling_frequency+dataset.initial_time;
z_axis= time*dataset.c0/2;
scan = linear_scan(scan.x_axis,z_axis);
u_k = zeros(484137,1);
v_k = zeros(484137,1);
X = zeros(1251,387);
lambd_k1 = zeros(484137,1);
lambd_k2 = zeros(1251,387);
objective = 1e6;
options = [];
options.display = 'none';
options.maxFunEvals = 100;
options.Method = 'lbfgs';
for i=1:maxiter
    %%%%%%%%%%%%%%%%%%%%%%%%%% update u %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % argmin_u  1/2*||b-u||^2 + beta1/2*||u - v + lambd_k1/beta1||^2 + beta2/2*||G - u + lambd_k2/beta2||^2
    %u_k = minFunc(@(u) FunCtion(b, Phi, u, v_k, lambd_k1, lambd_k2(:), beta1, beta2, G(:)), zeros(484137,1), options);
    u_k = (DAS_RF+beta2*(X(:)+lambd_k2(:)/beta2)+beta1*(v_k-lambd_k1/beta1))/(1+beta1+beta2);
    %%%%%%%%%%%%%%%%%%%%%%%%%% update G %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % argmin_G ||G||+ beta2/2*||G - u + lambd_k2/beta2||^2
    Gk = reshape(u_k,[1251 387])-lambd_k2/beta2;
    X = nuclear_norm_matrix(Gk,lambda/beta2,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%% update v %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % argmin_v  mu||v||_1 + beta1/2*||u - v + lambd_k1/beta1||^2
    %v_k = max(abs(u_k+lambd_k1/beta1)-mu/beta1,0).*sign(u_k+lambd_k1/beta1);
    % Replace the l1-norm soft thresholding with l0-norm hard thresholding
    v_k = max(abs(u_k+lambd_k1/beta1)-mu/beta1,0).*sign(u_k+lambd_k1/beta1);
    %%%%%%%%%%%%%%%%%%%%%%%%%% update lambd %%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambd_k1 = lambd_k1+beta1*(u_k-v_k);
    lambd_k2 = lambd_k2+beta2*(X-reshape(u_k,[1251 387]));
    %%%%%%%%%%%%%%%%%%%% stopping criterion %%%%%%%%%%%%%%%%%%%%%%%
    objective(i+1) = 0.5*norm(b-full(Phi*sparse(double(u_k))))^2 + mu * sum(abs(v_k))...
        +(beta1/2)*norm(u_k-v_k+(lambd_k1/beta1))^2 +...
        +(beta2/2)*norm(X(:)-u_k+(lambd_k2(:)/beta2))^2+lambda*sum(abs(X(:)));
    criterion(i) = abs(objective(i+1)-objective(i))/objective(i);
    disp(['iteration ',num2str(i),', criterion: ',num2str(criterion(i))])
    
    if ( criterion(i) < eps )
        u(:,i) = u_k;
        v(:,i) = v_k;
        lambd_k11(:,i) = lambd_k1;
        lambd_k22(:,i) = lambd_k2(:);        
        break
    end
    u(:,i) = u_k;
    v(:,i) = v_k;
    lambd_k11(:,i) = lambd_k1;
    lambd_k22(:,i) = lambd_k2(:);
end
objective = objective(2:end);
end

