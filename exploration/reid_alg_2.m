function [x,r] = reid_alg(A,b,tol)

% This is a really crappy implementation of Reid's algorithm
% for testing

n = size(A,1);
r = zeros(n,1);
x = zeros(n,1);

r = r + b;

while any(r >= 2*tol)
    j = find(r >= 2*tol, 1);
    x(j) = x(j) + (r(j) - (1.995)*tol);
    r = b - A*x;
end

