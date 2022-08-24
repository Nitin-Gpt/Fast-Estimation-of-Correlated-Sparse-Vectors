function [g_hat_corr,R_g,omega_inv_initial,x,omega_g,omega_inv_sum,sig,omega_c] = SBL_CORR_2(U,T,phi,y,noise_var, em_max_iter,em_thresh,D,M)
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
    if em_iter<10
        omega_g = (1/noise_var)* (phi') * phi + omega_c;
        x = inv(omega_g);
        omega_inv_initial = x;
        mu = (1/noise_var)*x*(phi')*y;
        [L,omega_inv_sum] = cal_L1(x);
        %L=round(3*D/4);
    else
        omega_g = (1/noise_var)* (phi') * phi + omega_c;
        [x,ind] = getInv(omega_g,x,D,L);
        omega_c_inv = inv(omega_c);
        %x = omega_c_inv - omega_c_inv * (phi') * inv(noise_var*eye(M)) * phi * omega_c_inv;
        mu = linsolve(omega_g,(1/noise_var)*(phi')*y);
        %mu = (1/noise_var)*x*(phi')*y;
    end
    
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
%ind;
sig = 1/noise_var;
omega_inv_sum = sort(ones(1,D)*abs(omega_inv_initial));
end