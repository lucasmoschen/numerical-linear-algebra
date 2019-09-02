function [B,P] = diag_dom(A)
    //função que tenta tornar uma matriz em uma matriz com diagonal estritamente dominante
    //B é a matriz com diagonal estritamente dominante e P é a matriz de permutações
    vetor_trocas = []
    for i=1:size(A,1)
        if A(i,i) <= norm(A(i,:),1) - abs(A(i,i))
            for j=1:size(A,2)
                if 
    end
endfunction
