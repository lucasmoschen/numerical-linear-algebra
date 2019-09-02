function [lambda,x,N] = Metodo_potencia_simetrico(A,x0,epsilon,M)
    //encontra o maior raio espectral de matrizes simétricas
    k = 1
    x = x0/norm(x0)
    y = A*x
    lambda = x'*y
    erro = norm((x - (y/norm(y,2))),2)
    x = (y/norm(y,2))
    k + 1
    while k <= M && erro > epsilon
        y = A*x
        lambda = x'*y
        if norm(y,2) == 0 then
            disp("A matriz tem autovalor 0")
        end
        erro = norm(abs(x)-abs(y/norm(y,2)),2)
        x = (y/norm(y,2))
        k = k + 1
    end
    N = k
    if k > M then
        disp("Número de iterações excedido")
    end
endfunction
