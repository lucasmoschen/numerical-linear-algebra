function singular_values(A,p)
% A is an image represeted by a matrix and p is vector with numbers between 0 and 1.
[U,S,V] = svd(double(A));

r = diag(S);
r = r(r > 0);

c = ceil(sqrt(size(p,1)));
f = floor(sqrt(size(p,1)));

for j = 1:size(p,1)
    s = max(1,floor(p(j)*size(r)));
    subplot(c, f, j), imshow(U(:,1:s)*S(1:s,1:s)*(V(:,1:s))', []);
end
end

