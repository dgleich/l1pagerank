%% Compute the examples for the graph from figure 1 in the Alg-Anti-Diff paper

%% Setup the graph
G = zeros(10);
G(1,2) = 1;
G(1,3) = 1;
G(2,3) = 1;
G(2,4) = 1;
G(2,6) = 1;
G(3,4) = 1;
G(3,5) = 1;
G(3,6) = 1;
G(3,9) = 1;
G(3,10) = 1;
G(4,5) = 1;
G(4,6) = 1;
G(5,6) = 1;
G(5,7) = 1;
G(6,7) = 1;
G(6,9) = 1;
G(6,10) = 1;
G(7,8) = 1;
G(8,9) = 1;
G(9,10) = 1;
G = G + triu(G)';
A = sparse(G);

S = 1:5;

%% Setup the cut problem
alpha = 0.5;

M = gmatrices(A); % extract matrices from a graph
n = M.n;

% Setup the vectors for the set S
Sv = zeros(M.n,1);
Sv(S) = 1;
Sbar = ~Sv;

%% Setup the localized cut graph
volS = sum(M.d.*Sv);

AS = [0              alpha*(M.D*Sv)'   0;
      alpha*(M.D*Sv)  A                alpha*(M.D*Sbar);
      0              alpha*(M.D*Sbar)' 0];

MS = gmatrices(AS);

%% Setup the PageRank problems
beta = 1/(1+alpha);
v = (M.d.*Sv)/volS;
tau = 1e-2;
kappa = tau*volS/beta;

%% Compare the PageRank solution to the 2-norm version

% PageRank
cvx_begin
    variable x(n+2);
    minimize norm(MS.Chalf*MS.Bu*x,2)
    subject to
        x(1) == 1;
        x(n+2) == 0;
cvx_end
z = x(2:end-1);
zn = z/volS;

pr = (alpha*M.D + M.K)\(alpha*v);
pr2 = ((eye(M.n) - beta*M.P')\((1-beta)*v))./M.d;

%% Show the solutions (should be all the same)
[zn pr pr2]

%% Print for the table given after our first theorem.
for i=1:10
    fprintf('%i ', M.d(i));
    fprintf('& %.4f ', pr2(i).*M.d(i));
    fprintf('& %.4f ',pr(i));
    fprintf('& %.4f ',zn(i)*volS);
    fprintf('\\\\ \n');
end

%% (See below for thie second table)

%% Solve the cut problem
cvx_begin
    variable x(n+2);
    minimize norm(MS.C*MS.Bu*x,1)
    subject to
        x(1) == 1;
        x(n+2) == 0;
cvx_end
zc = x(2:end-1);
zcn = zc/volS;

%% Show the PageRank vs. cut solution (should be different)
[zn pr zcn]

%% Solve the ACL problem
d2 = [0; M.d; 0];
cvx_begin
    cvx_precision high
    variable x(n+2);
    minimize 1/2*pow_pos(norm(MS.Chalf*MS.Bu*x,2),2) + kappa*sum(d2.*x)
    subject to
        x >= 0;
        x(1) == 1;
        x(n+2) == 0;
cvx_end
yacl = x(2:end-1);
xg = M.d.*yacl/volS;

%%
%[xr,rr] = reid_alg_pure_fixed(M.P, (1-beta)*v, M.d, beta, tau);
[xr,rr] = acl_method(M.P, (1-beta)*v, M.d, beta, tau, 1);

yr = (xr./M.d)*volS

%%
[xg xg - xr xr]

%% Print out the table
for i=1:10
    fprintf('%i ', M.d(i));
    fprintf('& %i ', round(zc(i)));
    fprintf('& %.4f ',pr(i)*volS);
    fprintf('& %.4f ',yacl(i));
    fprintf('\\\\ \n');
end

%% Print out the standard PageRank table
alpha1 = 3/17;
beta1 = 1/(1+alpha1);
svec1 = ones(n,1)/n;
sbar1 = M.d - svec1;
svol1 = n;
AS1 = [0              alpha*(svec1)'   0;
      alpha*(svec1)  A                alpha*(sbar1);
      0              alpha*(sbar1)' 0];

MS1 = gmatrices(AS1);

cvx_begin
    variable x(n+2);
    minimize norm(MS1.Chalf*MS1.Bu*x,2)
    subject to
        x(1) == 1;
        x(n+2) == 0;
cvx_end
z1 = x(2:end-1);
zn1 = z1;

prfull = ((eye(M.n) - beta*M.P')\((1-beta)*ones(n,1)/n))./M.d;

[zn1 prfull] % these should be equal

%% And now print the actual table
for i=1:10
    fprintf('%i ', M.d(i));
    fprintf('& %.4f ', prfull(i)*M.d(i));
    fprintf('& %.4f ', z1(i));
    fprintf('\\\\ \n');
end

%% Other stuff
% Reid

%y = x;
%xg = M.d.*y/volSbar;
%[xg xg - xr xr]

%% CVX variations on the QP
Q = alpha*M.D + M.K;
cvx_begin
    variable x(n);
    minimize 1/2*quad_form(x,Q) - alpha*sum(x.*M.d.*Sv) + kappa*sum(M.d.*x)
    subject to
        x >= 0;
cvx_end
yacl = x;
% minimize 1/2*pow_pos(norm(MS.Chalf*MS.Bu*x,2),2) + kappa*sum(MS.d.*x)
% minimize 1/2*quad_form(MS.Bu*x,MS.C) + kappa*sum(MS.d.*x)


