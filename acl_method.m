function [x,r,nsteps,edges] = acl_method(P,v,d,alpha,tol,rho)
% ACL_METHOD Solve PageRank via the ACL method
%
% [x,r] = acl_method(P,v,d,alpha,tol,rho)
%   finds x such that 0 <= v + alpha*P* x - x <= tol*d
% where each step involves increaing x(j) by r(j)-rho*d(j)*tol.
% thus, setting rho = 1 results in the sparsest vector possible
% and setting rho = 0 results in the densest vector possible. 
%
% When rho = 1, then the method also solves a particular 1-norm regularized
% problem.

% This is a really crappy implementation of Reid's algorithm
% for testing

n = size(P,1);
x = zeros(n,1);
r = v;

inQ = zeros(n,1);
Q = zeros(n,1);
qhead = 1;
qtail = 1;

for i = find(r)'
    if r(i) > d(i)*tol 
        inQ(i) = 1;
        Q(qtail) = i;
        qtail = mod(qtail,n) + 1;
    end
end

Pt = P';
iter = 1;
nedges = 0;

while qtail-qhead + n*(qtail < qhead) > 0
    j = Q(qhead);
    qhead = mod(qhead,n) + 1;
    inQ(j) = 0;
    
    rj = r(j);
    incr = r(j) - rho*d(j)*tol;
    assert(incr > 0);
    x(j) = x(j) + incr;

    r(j) = rho*d(j)*tol;
    for i = find(Pt(:,j))'
        assert(i ~= j);
        nedges = nedges + 1;
    	r(i) = r(i) + alpha*incr*Pt(i,j);
        if ~inQ(i) && r(i) > d(i)*tol 
            inQ(i) = 1;
            Q(qtail) = i;
            qtail = mod(qtail,n) + 1;
            assert(qtail ~= qhead); 
        end
    end
    iter = iter + 1;
    if mod(iter,n) == 0
        fprintf('step = %4i, sumr = %9.3e, qlen = %8i\n', ...
            iter, sum(r), qtail-qhead + n*(qtail < qhead));
    end
end

nsteps = iter;
edges = nedges;
