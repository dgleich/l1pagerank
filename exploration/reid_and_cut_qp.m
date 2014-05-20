%% Does Reid's algorith solve the cut QP?
%
% If so, how do the values of kappa and tau relate?
%
% x - solution of Reid/PageRank (sums to 1 ...)
% y - the soution to the cut QP (D^{-1} x)


%graph = 'four-clusters';
graph = 'dolphins';
%graph = 'netscience-cc';

[A,xy] = load_graph(graph);

%%
M = gmatrices(A);
n = M.n;

alpha = 0.85;
beta = 1/(1+alpha);

Sbar = zeros(M.n, 1);
Sbar(1) = 1;
S = ~Sbar;

AS = [0              alpha*(M.D*S)'   0;
      alpha*(M.D*S)  A                alpha*(M.D*Sbar);
      0              alpha*(M.D*Sbar)' 0];

  
MS = gmatrices(AS);

volSbar = sum(M.d.*Sbar);

tau = 1e-4;
kappa = tau*volSbar/beta;


%% Solve Reid's problem

[xr,rr] = reid_alg_pure(M.P, (1-beta)*Sbar, M.d, beta, tau, 1-1e-5);

%%
% This one uses rho = 1 and a GS algorithm
[xr,rr] = reid_alg_pure_fixed(M.P, (1-beta)*Sbar, M.d, beta, tau);

%% Here is our residual
rrr = (1-beta)*Sbar - xr + beta*M.P'*xr;

sr = tau*M.d - rrr; % this is positive
assert(all(sr >= 0))

%% Conver to y
y = xr./M.d;
r2 = (1-beta)*Sbar - (M.D - beta*M.A)*y;
s2 = tau*M.d - r2;
assert(all(s2 >= 0));

%% 
r3 = (1-beta)/beta*Sbar - (M.D/beta - M.A)*y;
s3 = tau/beta*M.d - r3;
assert(all(s3 >= 0));

%%
r4 = alpha*Sbar - (alpha*M.D + M.K)*y;
s4 = tau/beta*M.d - r4;
assert(all(s4 >= 0));

%%

(M.d.*tau - rr).*xr

%% Check out variations of this

(tau/beta)*M.d - alpha*Sbar + (alpha*M.D + M.K)*xr./M.d

%% Try out our proposed solution
z = volSbar*xr./M.d;
r = alpha*M.D*Sbar - (alpha*M.D + M.K)*z;
s = tau*volSbar/beta*M.d - r;

%% Figure out why the LP for Reid and the QP do the same thing.
Q = alpha*M.D + M.K;

cvx_begin
    variable x(n);
    minimize 1/2*quad_form(x,Q) - alpha*sum(x.*M.d.*Sbar) + kappa*sum(M.d.*x)
    subject to
        x >= 0;
cvx_end

% cvx_begin
%     variable x(n+2);
%     minimize norm(MS.Chalf*MS.Bu*x,2) + tau*sum(MS.d.*x)
%     subject to
%         x >= 0;
%         x(1) == 0;
%         x(n+2) == 1;
% cvx_end
% y = x(2:end-1);

y = x;
xg = M.d.*y/volSbar;
[xg xg - xr xr]

