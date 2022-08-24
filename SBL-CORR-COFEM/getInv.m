function [omega_inv,ind] = getInv(omega_g,omega_inv_old,D,L)
    omega_inv_sum = ones(1,D)*abs(omega_inv_old);
    [B,ind] = maxk(omega_inv_sum,L);
    A = zeros(D,L);
    for i = 1:L
        A(ind(i),i) = 1;
    end
    
    x = linsolve(omega_g,A);
    omega_inv = zeros(D,D);
    for i = 1:L
        omega_inv(:,ind(i)) = x(:,i);
    end
end