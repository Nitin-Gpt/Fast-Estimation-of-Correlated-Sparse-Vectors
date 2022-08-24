function [ R_R ] = exponentialcorrelation( N,l_R_R )
R_R = zeros(N,N);
for i=1:N
        for j=1:N
            R_R(i,j) = l_R_R^(abs(i-j));
        end
end
end