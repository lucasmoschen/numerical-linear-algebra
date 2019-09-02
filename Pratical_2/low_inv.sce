function Low_inv = inverseLmaisD(A)
    n = size(A,1)
    Low_inv = zeros(n,n)
    for j=1:n
        Low_inv(j,j) = 1/A(j,j)
        for i=j+1:n
            Low_inv(i,j) = (-1/A(i,i))*A(i,1:i-1)*Low_inv(1:i-1,j)
        end
    end
endfunction
