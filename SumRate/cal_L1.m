function [L,omega_inv_sum] = cal_L1(omega_inv)
    [D,D] = size(omega_inv);
    omega_inv_sum = sort(ones(1,D)*abs(omega_inv));
    %omega_sum = omega_inv_sum*ones(D,1);
    %omega_sum = omega_sum/D;
    L=0;
    for i=1:D
        if omega_inv_sum(i)>omega_inv_sum(D)*0.15
            L = L+1;
        end
    end
    if L>D
        L=D;
    end
end