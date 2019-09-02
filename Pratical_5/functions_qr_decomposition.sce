function [Q,R] = gram_schimidt(A)
    Q = zeros(A)
    R = zeros(size(Q,2),size(Q,2))
    for j = 1:size(A,2)
        vj = A(:,j)
        for i = 1:j-1
            R(i,j) = (A(:,j))'*Q(:,i)
            vj = vj - R(i,j)*Q(:,i)
        end
        R(j,j) = norm(vj)
        Q(:,j) = vj/R(j,j)
    end
endfunction

function [Q,R] = gram_schimidt_modificado(A)
    Q = zeros(A)
    R = zeros(size(Q,2),size(Q,2))
    for j = 1:size(A,2)
        vj = A(:,j)
        for i = 1:j-1
            R(i,j) = vj'*Q(:,i)
            vj = vj - R(i,j)*Q(:,i)
        end
        R(j,j) = norm(vj)
        Q(:,j) = vj/R(j,j)
    end
endfunction

function [U,R] = householder(A)
    R = A
    U = zeros(A)
    for i = 1:size(A,2)
       u = R(i:$,i)
       if u(1,1) > 0 then
        s = 1
       else
        s = -1
       end
       u(1,1) = u(1,1) + s*norm(u)
       U(i:$,i) = u/norm(u)
       H = eye(size(u,1),size(u,1)) - (2*u*u')/(u'*u)  // u ainda não é normado, pois 
       R(i:$,i:$) = H*R(i:$,i:$)                       //preferi não fazer a raiz quadrada
    end                              
endfunction

function x = resolve_sistema_escalona(R,c)
    C = [R, c]
    n = size(R,2)
    x=zeros(n,1);
    x(n)=C(n,n+1)/C(n,n);
    for i=n-1:-1:1
        x(i)=(C(i,n+1)-C(i,i:n)*x(i:n))/C(i,i);
    end
endfunction
