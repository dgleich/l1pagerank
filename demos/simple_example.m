%% Load the graph
load fig1ex.mat
M = gmatrices(A); % extract matrices from a graph, L, D, K, W
n = size(A,1);
Sv = zeros(M.n,1);
Sv(S) = 1;
Sbar = ~Sv;
volS = sum(M.d.*Sv);

%%
alpha = 0.5; beta = 1/(1+alpha);
tau = 1e-2; kappa = tau*volS/beta;
v = (M.d.*Sv)/volS;

%%
AS = [0              alpha*(M.D*Sv)'   0;
      alpha*(M.D*Sv)  A                alpha*(M.D*Sbar);
      0              alpha*(M.D*Sbar)' 0];
MS = gmatrices(AS);


%%
% PageRank <-> electrical flow
cvx_begin
    variable x(n+2);
    minimize norm(MS.Chalf*MS.Bu*x,2)
    subject to
        x(1) == 1;
        x(n+2) == 0;
cvx_end
z = x(2:end-1);
zn = z/volS;
% Compute PageRank
ztrue = ((eye(M.n) - beta*M.P')\((1-beta)*v));
[zn ztrue./M.d]

%%
% push method <-> 1 norm regularization
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
[xr,rr] = acl_method(M.P, (1-beta)*v, M.d, beta, tau, 1);
[xr,xg]
