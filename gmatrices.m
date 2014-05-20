function M = gmatrices(A,diagok)
% GMATRICES Compute all the matrices derived from a graph possible.
%
% M.n - number of vertices
% M.nedge - number of edges

if nargin == 1
    diagok = false;
end

M.A = A;

At = triu(A);

% basic stats
n = size(A,1);
nedge = nnz(At);
vol = nnz(A);

if ~diagok
    assert(all(diag(A) == 0)); % check that we have no diagonal entries
end

M.n = n;
M.nedge = nedge;
M.vol = vol;

% compute all the different degree matrices and vectors
     d = full(sum(A,2));             % degree vector
    dn = full(spfun(@(x) 1./x, d));  % inverse degree vector
 dhalf = sqrt(d);                    % square root of degrees
dnhalf = sqrt(dn);                   % inverse square root of degrees
     D = diag(sparse(d));
    Dn = diag(sparse(dn));
 Dhalf = diag(sparse(dhalf));
Dnhalf = diag(sparse(dnhalf));

M.d = d;
M.dn = dn;
M.dhalf = dhalf;
M.dnhalf = dnhalf;
M.D = D;
M.Dn = Dn;
M.Dhalf = Dhalf;
M.Dnhalf = Dnhalf;


% Compute the combinatorial matrices
[ei ej ev] = find(triu(A,0));

C = diag(sparse(ev));
Chalf = diag(sparse(sqrt(ev)));
c = full(ev);
chalf = full(sqrt(ev));

Bu = sparse([(1:nedge)'; (1:nedge)'], [ei ej], [ones(nedge,1); -1*ones(nedge,1)], nedge, n);
B = Chalf*Bu;
K = laplacian(A);

M.edges = [ei ej];

M.Bu = Bu; % unweighted incidence
M.B = B; % capacity weighted incidence
M.K = K; % combinatorial Laplacian
M.C = C; % edge weight matrices and vectors
M.Chalf = Chalf;
M.c = c;
M.chalf = chalf;


% Compute the normalized matrices
P = normout(A);
L = nlaplacian(A);
W = wadjacency(A);
N = B*Dnhalf;

M.P = P;
M.L = L;
M.W = W;
M.N = N;

