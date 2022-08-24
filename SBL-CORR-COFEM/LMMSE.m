function [h_hat] = LMMSE(R1,noise_var,W,y,M)

h_hat=R1*W'*(W*R1*W'+noise_var*eye(M))^(-1)*y;

end