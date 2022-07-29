function [c,kappa] = interpolate(x, f)
    n = size(x, 2);

    M = zeros(n,n);
    for j = 1:n
        for i = 1:n
            M(i,j) = cos((j-1)*acos(x(i)));
        end
    end

    b = zeros(n,1);
    for i = 1:n
        b(i,1) = f(x(i));
    end
    
    c = linsolve(M,b);
    kappa = cond(M);
end