function [x,r] = reid_alg_pure_fixed(P,v,d,alpha,tol)

% This is a really crappy implementation of Reid's algorithm
% for testing

n = size(P,1);
x = zeros(n,1);

r = v;

Pt = P';
iter = 1;

while any(r > d*tol)
    for j = 1:n
        if r(j) <= d(j)*tol, continue; end
        % What's actually done in Reid
        %x(j) = x(j) + (1-alpha)*(r(j) - rho*d(j)*tol);
        %r = v - x/(1-alpha) + alpha/(1-alpha)*P'*x
        % but this is the same

        incr = r(j) - d(j)*tol;
        assert(incr > 0);
        x(j) = x(j) + incr;

        r(j) = d(j)*tol;
        for i = find(Pt(:,j))
            r(i) = r(i) + alpha*incr*Pt(i,j);
        end
    end
    iter = iter + 1;
    fprintf('iter = %8i  sumr = %.16e\n', iter, sum(r));
    
    %r = v - x + alpha*P'*x;
    
end

