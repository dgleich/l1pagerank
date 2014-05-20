function [x,r] = reid_alg(A,b,tol,rho)

% This is a really crappy implementation of Reid's algorithm
% for testing

n = size(A,1);
r = zeros(n,1);
x = zeros(n,1);

r = r + b;

while any(r >= tol)
    j = find(r >= tol, 1);
    x(j) = x(j) + (r(j) - (rho)*tol);
    r = b - A*x;
end

