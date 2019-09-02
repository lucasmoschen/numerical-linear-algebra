


function [xk, dif, k, residuo] = Jacobi_solver(A,b,x0,E,M,norma)
    //A é a matriz
    //b é o vetor
    //x0 é aproximação inicial 
    //E é a tolerancia do erro
    //M é o número máximo de iterações
    //norma pode se 1, 2 ou %inf
    
    //xk é a solução aproximada.
    //dif é ||nk - nk-1||
    //k é o número de iterações
    //residuo é ||b - Axk||
    //criterio é o critério de parada
    n = size(A,1)
    if ~(n == size(A,2)) then
        error "A matriz deve ser quadrada"
    end

    D = diag(diag(A))
    D_inv = zeros(n,n)
    for i=1:n
        D_inv(i,i) = 1/D(i,i)
    end
    J = (-1)*D_inv*(A - D) //matriz do método de Jacobi L + U = A - D
    c = D_inv*b 
    if ~exists("norma", "local") then
        norma = 2
    end
    
    xk = J*x0 + c
    dif = norm(xk - x0, norma)
    residuo = norm(b - A*xk, norma)
    k = 1
    while dif >= E && k <= M
        xk_prev = xk
        xk = J*xk + c
        dif = norm(xk - xk_prev, norma)
        residuo = norm(b - A*xk, norma)
        k = k + 1
    end
    if dif < E then
        disp("||x_k - x_k-1|| < E: a norma da diferença entre os vetores é menor do que a tolerância")
        disp(" ")
    else
        disp("k>M: o número de iterações ultrapassou o máximo")
        disp(" ")
    end           

endfunction
