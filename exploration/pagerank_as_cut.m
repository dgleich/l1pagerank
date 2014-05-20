%% The relationships between PageRank and a min cut

graph = 'four-clusters';
%graph = 'dolphins';
%graph = 'netscience-cc';

[A,xy] = load_graph(graph);

M = gmatrices(A);
n = M.n;



%% When we convert a cut into PageRank, 
% the teleportation vector is like the opposite side of the cut. This
% is a matter of choice and definition.
% 
Sbar = zeros(M.n, 1);
Sbar(1) = 1;
S = ~Sbar;
alpha = 0.85;

AS = [0              alpha*(M.D*S)'   0;
      alpha*(M.D*S)  A                alpha*(M.D*Sbar);
      0              alpha*(M.D*Sbar)' 0];

MS = gmatrices(AS);

%beta = (1-alpha)/alpha;
beta = 1/(1+alpha);
Lpr = speye(n) - beta*M.W;

%% Check basic properties
norm(MS.L - MS.N'*MS.N,'fro')
norm(MS.L(2:end-1,2:end-1) - Lpr,'fro')


%% Setup the cut problem
cvx_begin
    variable x(n+2);
    minimize norm(MS.C*MS.Bu*x,1)
    subject to
        x >= 0;
        x(1) == 0;
        x(n+2) == 1;
cvx_end

%% Compute the max-flow
[flowval,cut,R,F] = max_flow(AS,1,n+2);
cut = (1-cut)/2;
norm(x-cut)

%% Setup the two-norm cut problem
cvx_begin
    variable x(n+2);
    minimize norm(MS.Chalf*MS.Bu*x,2)
    subject to
        x >= 0;
        x(1) == 0;
        x(n+2) == 1;
cvx_end

%% Setup the PageRank problem
xpr = (speye(n) - beta*M.P')\((1-beta)*M.D*Sbar/sum(M.D*Sbar));
xg = M.Dn*(xpr*sum(M.D*Sbar));
norm([0; xg; 1] - x)

%% Setup a cut for a PageRank problem with uniform teleportation
% The idea here is to use an in-flow such that the out-flow is
% equal to the teleportation vector and the sum of in + out is
% the degree vector.  

Sbar = zeros(M.n, 1);
Sbar(1) = 1;
Sbar(2) = 1;
S = ~Sbar;
alpha = 0.85;

S = M.d - Sbar;

AS = [0              alpha*(S)'   0;
      alpha*(S)      A                alpha*(Sbar);
      0              alpha*(Sbar)' 0];

MS = gmatrices(AS);
beta = 1/(1+alpha);

% Solve the two-norm cut
cvx_begin
    variable x(n+2);
    minimize norm(MS.Chalf*MS.Bu*x,2)
    subject to
        x >= 0;
        x(1) == 0;
        x(n+2) == 1;
cvx_end

% Compare to the PageRank solution
xpr = (speye(n) - beta*M.P')\((1-beta)*Sbar/sum(Sbar));
xg = M.Dn*(xpr*sum(Sbar));
norm([0; xg; 1] - x)

%% Check what happens if we restrict the norm of x vs. reid's alg
Sbar = zeros(M.n, 1);
Sbar(1) = 1;
S = ~Sbar;
alpha = 0.85;

AS = [0              alpha*(M.D*S)'   0;
      alpha*(M.D*S)  A                alpha*(M.D*Sbar);
      0              alpha*(M.D*Sbar)' 0];

MS = gmatrices(AS);

beta = 1/(1+alpha);

% Setup the two-norm cut problem with restricted norm of x
% the norm is 1.5 because x(n+2) = 1.
cvx_begin
    variable x(n+2);
    minimize norm(MS.Chalf*MS.Bu*x,2)
    subject to
        x >= 0;
        x(1) == 0;
        x(n+2) == 1;
        sum(x) == 1.7228; % four-clusters, tol=1e-2
cvx_end

%xpr = (speye(n) - beta*M.P')\((1-beta)*M.D*Sbar/sum(M.D*Sbar));
%xg = M.Dn*(xpr*sum(M.D*Sbar));
%norm([0; xg; 1] - x)

xpr = reid_alg_3(speye(n) - beta*M.P', ((1-beta)*M.D*Sbar/sum(M.D*Sbar)), 1e-2);
xg = M.Dn*(xpr*sum(M.D*Sbar));
zz = [[0; xg; 1] x]
sum(zz)
norm(zz(:,1)-zz(:,2))

%% What type of values of kappa do we need?
% My theory is now suggesting that kappa > beta.

Sbar = zeros(M.n, 1);
Sbar(1) = 1;
S = ~Sbar;
alpha = 0.85;

AS = [0              alpha*(M.D*S)'   0;
      alpha*(M.D*S)  A                alpha*(M.D*Sbar);
      0              alpha*(M.D*Sbar)' 0];

MS = gmatrices(AS);

beta = 1/(1+alpha);

kappa = 0.30;
cvx_begin
    variable x(n+2);
    minimize pow_pos(norm(MS.Chalf*MS.Bu*x,2),2) + kappa*sum(x)
    subject to
        x >= 0;
        x(1) == 0;
        x(n+2) == 1;
cvx_end

%xpr = (speye(n) - beta*M.P')\((1-beta)*M.D*Sbar/sum(M.D*Sbar));
%xg = M.Dn*(xpr*sum(M.D*Sbar));
%norm([0; xg; 1] - x)

xpr = reid_alg_3(speye(n) - beta*M.P', ((1-beta)*M.D*Sbar/sum(M.D*Sbar)), 1e-2);
xg = M.Dn*(xpr*sum(M.D*Sbar));
zz = [[0; xg; 1] x]
sum(zz)
norm(zz(:,1)-zz(:,2))

