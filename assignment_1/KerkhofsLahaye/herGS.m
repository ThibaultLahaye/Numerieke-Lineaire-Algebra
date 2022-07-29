function [Q,R] = herGS(A)

    % A a matrix.
    % Returns Q an ortho gonal matrix and R an
    % upper triangular matrix such that A = QR.

    [m,n] = size(A);
     
    Q = zeros(m,n);
    R = zeros(n);
    
    for j=1:n                               % Gram-Schmidt orthogonalization
        v = A(:,j);                                 % y begins as column j of A
        for re=1:2
            for i=1:j-1
                h = Q(:,i)' * v;              % modify A(:,j) to v for more accuracy
                v = v - Q(:,i)*h;                % subtract to projection (q^T*a)*a
                R(j,j) = R(j,j) + h;
            end
        end                                         % v is now perpendicular to all of q
        R(j,j) = (v' * v) ^ (1/2);
        Q(:,j) = v / R(j,j);                        % normalize v to be the next unit vector q
    end
    
%         y = A(:,j);                                 % y begins as column j of A
%         for i=1:j-1
%             R(i,j) = Q(:,i)' * A(:,j);              % modify A(:,j) to v for more accuracy
%             y = y - R(i,j) * Q(:,i);                % subtract to projection (q^T*a)*a
%         end                                         % v is now perpendicular to all of q
%         R(j,j) = (y' * y) ^ (1/2);
%         Q(:,j) = y / R(j,j); 

    return
end