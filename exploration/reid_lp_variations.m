%% In this script, we study variations on Reid's LP 
% to understand which constraints are required on the LP in order to solve
% it.  It's basically a guess and check type computational exercise to
% guess the right result about the necessary constraints.

addpath('~/Dropbox/matlab-extern/cvx/')
cvx_setup;

%%
%G = load_graph('dolphins');
G = load_graph('karate');
%G = load_graph('four-clusters');
%G = load_graph('netscience-cc');

%%
n = size(G,1);
L = nlaplacian(G);
Ls = speye(n) - L;
alpha = 0.85;
A = speye(n) - alpha*Ls;
d = sum(G,2);
dhalf = full(sqrt(d));
b = zeros(n,1); b(1) = (1-alpha)/dhalf(1);
v = zeros(n,1); v(1) = 1;
P = normout(G);

%%
% The winning formulation!
% We slightly parameterize Reid's algorithm, and add a factor of rho, to
% control how much residual we leave at each node.  Rho = 0 leaves no extra
% residual at each node.  Rho = 1 leaves tol residual at each node. The LP
% solves the case with rho = 1.  In general, it's better to use rho =
% 0.9999.
tol = 3e-3;
cvx_begin
  variables y(n)
  minimize 2*norm(y,1)
  subject to 
    norm(A*y-b,inf) <= tol;
    y >= 0;
cvx_end
y = full(y);
r = b - A*y;
x = dhalf.*y;
sum(x)
min(y)
full([x(1:10) r(1:10) d(1:10)])
[sum(x > 1e-8) sum(r > 1e-8) ]
[yr,rr] = reid_alg_4(A,b,tol,1-1e-5); xr = dhalf.*yr;
[sum(xr > 1e-8) sum(rr > 1e-8)]
plot(x,xr,'.');
norm(x-xr,'inf')

%%
% We need the constraints 
% y >= 0, otherwise the Karate club network fails.


%% Find better variations of Reid's LP
cvx_begin
    variables y(n)
    minimize sum(y)
    subject to 
        b - A*y <= tol;  % required, otherwise we don't have our upper bound
        y >= 0;          % required, otherwise, solution is unbounded below
cvx_end    
y = full(y);
r = b - A*y;
x = dhalf.*y;
norm(x-xr,'inf')

%%
% We don't seem to need the constraint b-A*y >= 0, why not?
% This should work with the other variation of x too...

cvx_begin
    variables x(n)
    minimize sum(dhalf.*x)
    subject to 
        (1-alpha)*v + alpha*P'*x - x <= dhalf.*tol;  % required, otherwise we don't have our upper bound
        x >= 0;          % required, otherwise, solution is unbounded below
cvx_end    
norm(x-xr,'inf')

%%
% Some quick analysis showed that we may need this, if the diagonal of G
% has any non-zeros.  So I'm going to make them all non-zero!
%
% UPDATE: This experiment has been eliminated as we can handle it with the
% result below. This result shows that the residual should always be
% non-negative.
%
% G = G - diag(diag(G));
% G = G + speye(size(G));
% n = size(G,1);
% L = nlaplacian(G);
% Ls = speye(n) - L;
% alpha = 0.85;
% A = speye(n) - alpha*Ls;
% d = sum(G,2);
% dhalf = full(sqrt(d));
% b = zeros(n,1); b(1) = (1-alpha)/dhalf(1);
% v = zeros(n,1); v(1) = 1;
% P = normout(G);
% 
% cvx_begin
%     variables x(n)
%     minimize sum(dhalf.*x)
%     subject to 
%         (1-alpha)*v + alpha*P'*x - x <= dhalf.*tol;  % required, otherwise we don't have our upper bound
%         x >= 0;          % required, otherwise, solution is unbounded below
% cvx_end    
% [yr,rr] = reid_alg_4(A,b,tol,1-1e-5); xr = dhalf.*yr;
% 
% norm(x-xr,'inf')

%%
% Okay, turns out that also wasn't the case (see update), the lemma is the following
% Suppose that any component of the residual b - A*x is smaller than zero.
% then we know that b(i) - [(I - alpha W)*x](i) < 0 for some i.
% For this i, we have x(i) > b(i) + alpha w_i^T*x >= 0.
% So if we set x(i) to (1-alpha w_ii)^{-1} b(i) + alpha w_i^T*x, 
% then we will reduce the objective and ...
%    b - (I - alpha W)(x - z e_i) 
%         = b - (I - alpha W)*x + z e_i - alpha z W e_i
%         <= tol*e - alpha z W e_i + z e_i
% So for all components except i, we have:
%    b - (I - alpha W)(x - z e_i) <= tol*e - alpha z W e_i
% which will satisfy the constaint.
% In the ith component, we have:
%    b(i) - x(i) + alpha w_i^T*x + z - alpha z w_ii
% but the term z = x(i) - b(i) - alpha w_i^T*x
% which will make this term zero, and satisfy the constraint.
% So this is redundant.

%% Get all the constants right.
% Now we want to compare against the pure version of Reid's algorithm in
% order to make sure we get each constant correct.

[xr,rr] = reid_alg_pure(P,v,d,alpha,tol,1-1e-5);
cvx_begin
    variables x(n)
    minimize sum(x)
    subject to 
        v + alpha*P'*x - x <= d.*tol;  % required, otherwise we don't have our upper bound
        x >= 0;          % required, otherwise, solution is unbounded below
cvx_end    
x = (1-alpha)*x;
norm(x-xr,'inf')

%% Laplacian version
%beta = 1/(1-alpha);
beta = (1-alpha)/alpha;
[xr,rr] = reid_alg_pure(P,v,d,alpha,tol,1-1e-5);
cvx_begin
    variables z(n)
    minimize sum(d.*z)
    subject to 
        beta*v - (beta*diag(d) + (diag(d) - G))*z <= beta*d.*tol;  % required, otherwise we don't have our upper bound
        z >= 0;          % required, otherwise, solution is unbounded below
cvx_end    
x = full(d).*z;
norm(x-xr,'inf')

%% Equivalence between Laplacian and cut LP
% In a solution of the LP, why are these two terms equal?
% (I derived this based on the cut LP...
%   it's basically saying the solution satisfies complementary slackness)

beta*z'*v - z'*(beta*diag(d) + (diag(d) - G))*z
beta*tol*d'*z


%% Does the objective vector matter?

beta = (1-alpha)/alpha;
[xr,rr] = reid_alg_pure(P,v,d,alpha,tol,1-1e-5);
rv = rand(n,1);
cvx_begin
    variables z(n)
    minimize sum(rv.*z)
    subject to 
        beta*v - (beta*diag(d) + (diag(d) - G))*z <= beta*d.*tol;  % required, otherwise we don't have our upper bound
        z >= 0;          % required, otherwise, solution is unbounded below
cvx_end    
x = full(d).*z;
norm(x-xr,'inf')

%% Seems not.

%% What about negative entries?

beta = (1-alpha)/alpha;
[xr,rr] = reid_alg_pure(P,v,d,alpha,tol,1-1e-5);
rv = rand(n,1);
cvx_begin
    variables z(n)
    minimize sum(rv.*z)
    subject to 
        beta*v - (beta*diag(d) + (diag(d) - G))*z <= beta*d.*tol;  % required, otherwise we don't have our upper bound
        z >= 0;          % required, otherwise, solution is unbounded below
cvx_end    
x = full(d).*z;
norm(x-xr,'inf')