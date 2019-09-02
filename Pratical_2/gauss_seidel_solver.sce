//Por razões de exatidão no que desejo fazer, calculo a inversa de L + D manualmente. inv(L+D) muitas vezes retornou valores com um erro, mesmo que pequeno. 
function Low_inv = inversaLmaisD(A)
    n = size(A,1)
    Low_inv = zeros(n,n)
    for j=1:n
        Low_inv(j,j) = 1/A(j,j)
        for i=j+1:n
            Low_inv(i,j) = (-1/A(i,i))*A(i,1:i-1)*Low_inv(1:i-1,j)
        end
    end
endfunction

function [xk, dif, k, residuo] = GaussSeidel_solver(A,b,x0,E,M,norma)
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
    
    //calculo de (L + D)^(-1)
    Low_inv = inversaLmaisD(A)
    U = triu(A,1)
    GS = (-1)*Low_inv*U
    c = Low_inv*b
       
    if ~exists("norma", "local") then
        norma = 2
    end
    
    xk = GS*x0 + c
    dif = norm(xk - x0, norma)
    residuo = norm(b - A*xk, norma)
    k = 1
    while dif >= E && k <= M
        xk_prev = xk
        xk = GS*xk + c
        dif = norm(xk - xk_prev, norma)
        residuo = norm(b - A*xk, norma)
        k = k + 1
    end
    if dif < E then
        disp("||x_k - x_k-1|| < E: a norma da diferençça entre os vetores é menor do que a tolerância")
        disp(" ")
    else
        disp("k>M: o número de iterações ultrapassou o máximo")
        disp(" ")
    end           

endfunction
