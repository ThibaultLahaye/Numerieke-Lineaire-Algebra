function [v] = evalCheb(a,x)
    n = size(a,2)-1;
    
    lambda = zeros(1,n+2);
    lambda(:, 1:2) = 1;
    lambda(1,3:end) = 2;
    
    v = zeros(1,size(x,2));
    for i = 1 : size(x,2)
        b = zeros(1, n+3);
        for k = n : -1 : 0
            b(k+1) = a(k+1) + lambda(k+2)*x(i)*b(k+2) - b(k+3);
        end
        v(i) = b(1);
    end
end