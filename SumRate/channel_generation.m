function [h,Rh] = channel_generation(D,T,cov_theta,A,N)

% ind_g=randperm(D,L);
% gamma=zeros(D,1);
% for i=1:L
%     gamma(ind_g(i))=1;
% end
% Gamma=diag(gamma);
% cov_theta=sqrt(Gamma)*U*sqrt(Gamma);
% G=(mvnrnd(zeros(D,1),0.5*cov_theta,T)+1i*mvnrnd(zeros(D,1),0.5*cov_theta,T));
G1 = sqrt(1/2)*(randn(N,T) + sqrt(-1)*randn(N,T));
Rh=A*cov_theta*A';
h=sqrtm(Rh)*G1;
% A=dftmtx(N)/sqrt(N);
% A1=dftmtx(D)/sqrt(D);
% x=randperm(D,N);
% A=A1(x,:);
% h=A*G.';
% disp(g);
end