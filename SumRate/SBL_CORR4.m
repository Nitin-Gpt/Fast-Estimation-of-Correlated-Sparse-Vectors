function [g_hat_corr,R_g] = SBL_CORR4(U,T,phi,y,noise_var, em_max_iter,D,M,W)

em_thresh=0.001;
        em_iter = 1;
        condition = 1;
        %formatSpec = 'iteration: %d, error: %f';
        c = ones(D, 1); 
        Uinv = (U)^(-1);
        omega_c = diag(c)*Uinv*diag(c);
        cov_c = inv(omega_c);
        mu_g = ones(D,T);
        tau_p=10;
       % phitY = phi'*y;
        while(em_iter <= em_max_iter && condition == 1)
            mu =  mu_g;
%             R_g=zeros(D,D);
%             for r=1:T
%             %E-step
             cov_g = cov_c - cov_c*phi'*((noise_var*eye(M) + phi*cov_c*phi')\(phi*cov_c));
            mu_g = cov_g*phi'*y/noise_var;
%             R_g = R_g + (cov_g +(mu_g(:,r)*mu_g(:,r)'))/T;
%             end 
            %cov_g = cov_c - (tau_p)*cov_c*phi'*((noise_var*(W*W') + (tau_p)*phi*cov_c*phi')\(phi*cov_c));
            %mu_g = cov_g*phi'*(noise_var*(W*W'))^(-1)*y;
            R_g = cov_g +(mu_g*mu_g')/T;
            %M-step
            c = ((real(Uinv.*(R_g.'))))\(1./c);
            cov_c = diag(1./c)*U*diag(1./c);
            x=mu_g-mu;
%             x=x(:);

            % condition
            condition = (norm(x,'fro')) > em_thresh;
            em_iter = em_iter+1;
            %disp(em_iter);
        end
        
        g_hat_corr=mu_g;
end