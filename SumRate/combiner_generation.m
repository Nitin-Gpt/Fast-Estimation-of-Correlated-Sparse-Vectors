  function [Wr,Hr,Hr_hat] = combiner_generation(R_avg,K,N,h_hat,M,T)
 %function [Wr] = combiner_generation(R_avg,K,N,h_hat,M,T,Hr_tilde)
Wr = zeros(K,M);
SIG = zeros(M,N);
R_avg = R_avg/(K); 
[V1,D1] = eig(R_avg);
  d = diag(D1);
  if(d(1)>=d(end))
      V = V1;
  else
      V = flip(V1,2);
  end   
 for i = 1:M
      SIG(i,i) = 1;
end 

Hr1 = SIG*V';
Hr2 = angle(Hr1);
Hr = (1/sqrt(N))*exp(sqrt(-1)*Hr2);   
Hr_tilde =(Hr*Hr')\Hr;
% Hr =  eye(N);
% Hr_tilde = eye(N);
%% Digital precoding matrix
for k=1:K
Wr(k,:) = (Hr_tilde*h_hat(:,(k-1)*T+1))';
end
 Hr_hat=Hr_tilde'*Hr;
end