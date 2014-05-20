%% The relationships between PageRank and a min cut

%graph = 'four-clusters';
graph = 'dolphins';
%graph = 'netscience-cc';

[A,xy] = load_graph(graph);

M = gmatrices(A);
n = M.n;

%% Figure out why the LP for Reid and the QP do the same thing.
Sbar = zeros(M.n, 1);
Sbar(1) = 1;
S = ~Sbar;
alpha = 0.85;

AS = [0              alpha*(M.D*S)'   0;
      alpha*(M.D*S)  A                alpha*(M.D*Sbar);
      0              alpha*(M.D*Sbar)' 0];

MS = gmatrices(AS);

beta = 1/(1+alpha);

kappa = 0.01;

%Q = (MS.Chalf*MS.Bu)'*MS.Chalf*MS.Bu;
Q = alpha*M.D + M.K;

cvx_begin
    variable x(n+2);
    minimize norm(MS.Chalf*MS.Bu*x,2) + kappa*sum(MS.d.*x)
    %minimize quad_form(x,Q) + kappa*sum(MS.d.*x)
    subject to
        x >= 0;
        x(1) == 0;
        x(n+2) == 1;
cvx_end

xg = x(2:end-1);

% compute s, the lack in the constraint
s = 2*(alpha*(M.D) + (M.L))*xg - 2*alpha*(M.D*Sbar) + kappa*M.d;

tol = max(s)


%xg = M.Dn*(xpr*sum(M.D*Sbar));
zz = [[0; xg; 1] x]
sum(zz)
norm(zz(:,1)-zz(:,2))



%%

Sbar = zeros(M.n, 1);
Sbar(1) = 1;
S = ~Sbar;
alpha = 0.85;

AS = [0              alpha*(M.D*S)'   0;
      alpha*(M.D*S)  A                alpha*(M.D*Sbar);
      0              alpha*(M.D*Sbar)' 0];

MS = gmatrices(AS);

beta = 1/(1+alpha);

kappa = 0.01;

Q = alpha*M.D + M.K;

cvx_begin
    variable x(n);
    minimize quad_form(x,Q) - 2*alpha*sum(x.*M.d.*Sbar) + kappa*sum(M.d.*x)
    subject to
        x >= 0;
cvx_end

xg = x;

xpr = M.D*(xg)/sum(M.D*Sbar);
r = (1-beta)*Sbar - xpr + beta*M.P'*xpr;

[xr,rr] = reid_alg_pure(M.P, Sbar, M.d, beta, max(r./M.d)/(1-beta), 1-1e-5);

%% Check the equivalence of objectives.
kappa/2*sum(x.*M.d) - alpha*sum(M.d.*Sbar.*x) - cvx_optval

%%
% Okay, so these are still equivalent, why can't I show this?
% for the QP, we should have a positive set of Lagrange multipiers.

s = 2*Q*x - 2*alpha*M.d.*Sbar + kappa*M.d

% So we do.

% And we should have complementary slackness
x'*s

% Which we do.
r = alpha*M.d.*Sbar - Q*x - kappa/2*M.d

%% 
% The reason I couldn't show this is because of a normalization
% difference with one of the codes. After that was fixed, things
% came out fine.

%%
% In theory, the solution of this LP should be the same
cvx_begin
    variable y(n);
    minimize kappa/2*sum(y.*M.d) - alpha*sum(M.d.*Sbar.*y)
    subject to
        2*Q*y - 2*alpha*M.d.*Sbar + kappa*M.d >= 0 ;
        y >= 0;
cvx_end
% Which is unbounded, of course.
