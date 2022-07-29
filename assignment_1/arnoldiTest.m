function [Q, H] = arnoldiTest(A, v, m)
    n = size(A, 1);
    H = zeros(m+1, m);
    Q = zeros(n, m+1);
    Q(:, 1) = v/norm(v);
    for k = 1:m
        v = A*Q(:, k);
        for j = 1:k
            H(j, k) = Q(:, j)'*v;
            v = v - H(j, k)*Q(:, j);
        end
        H(k+1, k) = norm(v);
        Q(:, k+1) = v/H(k+1, k);
    end
end