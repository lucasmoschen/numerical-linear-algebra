function [x, C, P]=Gaussian_Elimination_4(A, b)
    C=[A,b];
    n=size(C,1); //n é o número de linhas de C
    P = eye(n,n);
    for j=1:(n-1) //O pivô está na posição (j,j)
        [maior_pivo,indice] = max(abs(C([j:n],j)));
        indice = j + indice - 1
        if maior_pivo == 0
            return "O sistema é indeterminado"
        end
        C([j,indice],:) = C([indice,j],:); //troca as linhas j e a encontrada
        P([j,indice],:) = P([indice,j],:); //troca as linhas na matriz de permutação
           
        for i=(j+1):n
            C(i,j)=C(i,j)/C(j,j);
            C(i,j+1:n+1)=C(i,j+1:n+1)-C(i,j)*C(j,j+1:n+1);
        end
    end

    x=zeros(n,1);
    x(n)=C(n,n+1)/C(n,n);
    for i=n-1:-1:1
        x(i)=(C(i,n+1)-C(i,i:n)*x(i:n))/C(i,i);
    end
    C=C(1:n,1:n);
endfunction

function [lambda,x,N] = Metodo_potencia_inversa(A,x0,epsilon,alfa,M)
    if alfa == %inf then
        alfa = ((x0'*A)*x0)/(x0'*x0)
    end
    k = 1
    p = find(max(abs(x0))==abs(x0),1)
    x = x0/x0(p)
    Matrix = A - alfa*eye(A)               //Calculo da matriz (A - alfaI)
    [y] = Gaussian_Elimination_4(Matrix, x)
    p = find(max(abs(y))==abs(y),1)
    lambda = y(p)
    erro = norm((x - y/y(p)),%inf)
    x = y/y(p)
    k = k + 1 
    while k <= M && erro > epsilon
        y = Gaussian_Elimination_4(Matrix, -x)
        p = find(max(abs(y))==abs(y))
        lambda = y(p)
        if norm(y,2) == 0 then
            disp("A matriz tem autovalor 0")
        end
        erro = norm((x - y/y(p)),%inf)
        x = y/y(p)
        k = k + 1 
    end
    lambda = (1/lambda)+alfa
    N = k
    if k > M then
        disp("Número máximo de iterações excedido")
    end
    
endfunction

