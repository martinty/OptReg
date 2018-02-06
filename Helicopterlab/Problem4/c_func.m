function [c,ceq] = c_func(Z)
    alfa = 0.2; beta = 20; lambdat = 2*pi/3;
    ceq = [];
    c = zeros(length(Z)/8);

    for i = 1:length(C)
        xk = Z(6*(i-1)+1:6*i);
        c(i) = alfa * exp (-beta*(xk(1)-lambdat)^2)-xk(5);
    end 

end

