function [lambda,x,N] = Metodo_potencia(A,x0,epsilon,M)
    //determina uma aproximação do autovalor dominante da matriz A
    //supomos que os autovetores associados aos n autovetore são linearmente independentes.
    //essa suposição garante o sucesso, mas não é condição necessária.
    //supomos que |lambda1| > |lambda i|, para todo i.
    k = 1
    p = find(max(abs(x0))==abs(x0),1)
    x = x0/x0(p)
    y = A*x
    p = find(max(abs(y))==abs(y))
    lambda = y(p)
    erro = norm((x - y/y(p)),%inf)
    x = y/y(p)
    k = k + 1                         //já faço uma iteração anterior ao while para obter o primeiro valor de lambda e erro
    while k<=M && erro > epsilon
        y = A*x
        p = find(max(abs(y))==abs(y),1)
        lambda = y(p)
        if norm(y,2) == 0 then
            disp("A matriz tem autovalor 0")
        end
        erro = norm((x - y/y(p)),%inf)
        x = y/y(p)
        k = k + 1
    end
    N = k
    if k > M then
        disp("Número máximo de iterações excedido")
    end
endfunction
