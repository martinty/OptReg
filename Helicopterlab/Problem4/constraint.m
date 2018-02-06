function [C, Ceq] = constraint(z)

global N mx beta alpha lambda_t
C_temp = -Inf(N, 1);
lambda_k = z(1:mx:N*mx);
ek = z(5:mx:N*mx);

for n = 1:N
    C_temp(n) = alpha*exp(-beta*(lambda_k(n)-lambda_t)^2)-ek(n);
end

Ceq = [];
C = C_temp;

end

