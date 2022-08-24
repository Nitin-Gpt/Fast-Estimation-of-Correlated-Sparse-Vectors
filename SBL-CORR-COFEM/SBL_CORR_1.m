function [g_hat_corr,R_g,x,omega_g,omega_c,sigm] = SBL_CORR_1(U,T,phi,y,noise_var, em_max_iter,em_thresh,D,M,L,cov_theta)
%initialization
em_iter = 0;
condition = 1;
mu_old = 100*ones(D,T);
mu = ones(D,T);
c = ones(D, 1); 
U_inv = (U)^(-1);
omega_c = diag(c)*U_inv*diag(c);
while (em_iter <= em_max_iter && condition == 1)
    R_g = zeros(D,D);
    mu_old = mu;
    omega_g = (1/noise_var)* (phi') * phi + omega_c;
    x = linsolve(omega_g,eye(D));
    %omega_c_inv = inv(omega_c);
    %x = omega_c_inv - omega_c_inv * (phi') * inv(noise_var*eye(M) + phi*omega_c_inv*(phi')) * phi * omega_c_inv;
    %x = omega_c_inv - omega_c_inv * (phi') * inv(noise_var*eye(M)) * phi * omega_c_inv;
    mu = (1/noise_var)*x*(phi')*y;
    for r = 1:T
        R_g = R_g + (mu(:,r)*(mu(:,r)'))/T;
    end
    % M step
    R_g = R_g + x;

    c = linsolve(real(U_inv.*(R_g.')),(1./c));
    omega_c = diag(c)*U_inv*diag(c);

    % condition
    condition = norm(mu-mu_old,'fro') > em_thresh;

    em_iter = em_iter+1;
end
g_hat_corr = mu;
sigm = 1/noise_var;
end