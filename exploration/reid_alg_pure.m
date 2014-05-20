function [x,r] = reid_alg_pure(P,v,d,alpha,tol,rho)

% This is a really crappy implementation of Reid's algorithm
% for testing

n = size(P,1);
x = zeros(n,1);

r = v;

while any(r > d*tol)
    j = find(r > d*tol, 1);
    % What's actually done in Reid
    %x(j) = x(j) + (1-alpha)*(r(j) - rho*d(j)*tol);
    %r = v - x/(1-alpha) + alpha/(1-alpha)*P'*x
    % but this is the same
    x(j) = x(j) + (r(j) - rho*d(j)*tol);
    r = v - x + alpha*P'*x;
end

