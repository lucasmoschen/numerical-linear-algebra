//////////////////////////////////////////////////////////////////////////
//Variáveis de saída:
//x: solução do sistema Ax=b (assumimos que tal solução existe).
//C: Seja PA=LU a decomposição LU de PA.
//Então C(i,j)=L(i,j) para i>j e C(i,j)=U(i,j) para j>=i.
//////////////////////////////////////////////////////////////////////////
function X=Resolve_com_LU(C,B,P) //C é matriz da decomposição LU_Simples, P é a matriz de troca de linhas do Pivoteamento Parcial e B a matriz de vetores b.
    //Vou resolver PAX = LUX = PB -> LY = PB -> UX = Y
    if size(C,1) ~= size(C,2) then
        error("A matriz C deve ser quadrada (LU)")
    end
    if size(C,1) ~= size(B,1) then
        error("A matriz B e a matriz C devem ter o mesmo número de linhas")
    end
    n = size(B,1)
    m = size(B,2)
    if ~exists("P", "local") then
        P = eye(n,n)
    end
    PB = P*B
    Y = zeros(n,m)
    for i=1:n      //para cada linha da matriz L (parte inferior da C)
        Y(i,:) = PB(i,:)
        if i > 1 then
            Y(i,:) = Y(i,:) - C(i,1:i-1)*Y(1:i-1,:)
        end
    end
    
    X = zeros(n,m)
    for i=n:-1:1      //para cada linha da matriz U (parte superior da C)
        X(i,:) = Y(i,:)
        if i < n then
            X(i,:) = X(i,:) - C(i,i+1:n)*X(i+1:n,:)
        end
        X(i,:) = X(i,:)/C(i,i)
    end
endfunction
