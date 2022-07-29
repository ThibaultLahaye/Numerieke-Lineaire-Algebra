% function [Q,R] = modGS(A)
%     % QR factorization with MGS
%     
%     [m, n] = size(A);
%     Q = zeros(m,n);
%     R = zeros(n,n);
%     
%     for i=1:n
%         R(i,i) = norm(A(:,i));
%         Q(:,i) = A(:,i)/R(i,i);
%         R(i,i+1:n) = Q(:,i)'*A(:,i+1:n);
%         A(:,i+1:n) = A(:,i+1:n) - Q(:,i)*R(i,i+1:n);
%     end
% end

function [Q,R] = modGS(A)

    Q = zeros(size(A));
    R = zeros(size(A,2));

    for i = 1:size(A,2)

        R(i,i) = norm(A(:,i));
        Q(:,i) = A(:,i)/R(i,i);

        for j = i+1:size(A,2)
            R(i,j) = Q(:,i)'*A(:,j);
            A(:,j) = A(:,j) - R(i,j)*Q(:,i);
        end
    end

end

% function [Q,R]=modGS(A)
%     [m,n]=size(A);
%     V=A;
%     Q=zeros(m,n);
%     R=zeros(n,n);
%     for i=1:n
%         R(i,i)=norm(V(:,i));
%         Q(:,i)=V(:,i)/R(i,i);
%         for j=i+1:n
%             R(i,j)=Q(:,i)'*V(:,j);
%             V(:,j)=V(:,j)-R(i,j)*Q(:,i);
%         end
%     end
% end