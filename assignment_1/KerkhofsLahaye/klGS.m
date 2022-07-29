
function [Q,R] = klGS(A)

    % A a matrix.
    % Returns Q an ortho gonal matrix and R an
    % upper triangular matrix such that A = QR.

    Q = zeros(size(A));
    R = zeros(size(A,2));
    
    for j=1:size(A,2)                               % Gram-Schmidt orthogonalization
        y = A(:,j);                                 % y begins as column j of A
        for i=1:j-1
            R(i,j) = Q(:,i)' * A(:,j);              % modify A(:,j) to v for more accuracy
            y = y - R(i,j) * Q(:,i);                % subtract to projection (q^T*a)*a
        end                                         % v is now perpendicular to all of q
        R(j,j) = (y' * y) ^ (1/2);
        Q(:,j) = y / R(j,j);                        % normalize v to be the next unit vector q
    end
    
    return
end