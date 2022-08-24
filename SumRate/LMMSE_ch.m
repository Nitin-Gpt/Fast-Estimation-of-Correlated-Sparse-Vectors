function [h_hat] = LMMSE_ch(R1,noise_var,W,y,M)

% R1=array_response*cov_theta*array_response';
tau_p = 10;
%h_hat=R1*W'*(((tau_p)*W*R1*W')+(noise_var*(W*W')*eye(M)))^(-1)*y;
h_hat=R1*W'*(W*R1*W'+noise_var*eye(M))^(-1)*y;
end