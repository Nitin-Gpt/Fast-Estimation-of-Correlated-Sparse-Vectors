function [h,G, ind_g, array_response,cov_theta] = channel_generation(D,N,U,L,T)

ind_g=randperm(D,L);
gamma=zeros(D,1);
for i=1:L
    gamma(ind_g(i))=1;
end
Gamma=diag(gamma);
cov_theta=sqrt(Gamma)*U*sqrt(Gamma);
G=zeros(T,D);
for t=1:T
    g=(mvnrnd(zeros(D,1),0.5*cov_theta)+1i*mvnrnd(zeros(D,1),0.5*cov_theta));
    for i=1:D
        G(t,i)=g(i);
    end
end
% array_response=dftmtx(N)/sqrt(N);
A1=dftmtx(D)/sqrt(D);
x=randperm(D,N);
array_response=A1(x,:);
h=array_response*G.';
% disp(g);
end