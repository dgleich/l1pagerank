%% Show some examples of how the various parameter affect the netscience prob

% Initital generation
% graph = 'netscience-cc';
% A = load_graph(graph);
% xy = igraph_draw(A,'kk');
% save 'example_netscience_data.mat' A xy;

%%
load 'example_netscience_data'

M = gmatrices(A); % extract matrices from a graph
n = M.n;

%% Test a simple example
ei = zeros(n,1); ei(1) = 1;
x = (speye(n) - 0.85*M.P') \ ei;
example_netscience_fig(A,xy,max(80*x,12),sqrt(x))

%% Figure parameters
ns = 40; % node scale
nm = 6; % node min


%% Setup the parameters and cut graph
S = find(xy(:,1)>-2 & xy(:,1) <-1 & xy(:,2) < 0);
% Setup the vectors for the set S
Sv = zeros(M.n,1);
Sv(S) = 1;
Sbar = ~Sv;

volS = sum(M.d.*Sv);
alpha = 0.5;
beta = 1/(1+alpha);
tau = 1e-3;
kappa = tau*volS/beta;

AS = [0              alpha*(M.D*Sv)'   0;
      alpha*(M.D*Sv)  A                alpha*(M.D*Sbar);
      0              alpha*(M.D*Sbar)' 0];

MS = gmatrices(AS);
example_netscience_fig(A,xy,max(ns*Sv,nm),sqrt(Sv))
title(sprintf('%i nonzeros',nnz(Sv > 1e-5)));
xlim([-3.5432    1.8662]);
ylim([ -5.6786    2.5062]);
print(gcf,'netscience-figure-set.eps','-depsc2','-painters')
%% Compute PageRank
cvx_begin
    variable x(n+2);
    minimize norm(MS.Chalf*MS.Bu*x,2)
    subject to
        x(1) == 1;
        x(n+2) == 0;
cvx_end
z = x(2:end-1);
example_netscience_fig(A,xy,max(ns*z,nm),sqrt(z))
title(sprintf('%i nonzeros',nnz(z > 1e-5)));
xlim([-3.5432    1.8662]);
ylim([ -5.6786    2.5062]);
print(gcf,'netscience-figure-pr.eps','-depsc2','-painters')
%% Compute cut
cvx_begin
    variable x(n+2);
    minimize norm(MS.C*MS.Bu*x,1)
    subject to
        x(1) == 1;
        x(n+2) == 0;
cvx_end
z = x(2:end-1);
example_netscience_fig(A,xy,max(ns*z,nm),sqrt(z))
title(sprintf('%i nonzeros',nnz(z > 1e-5)));
xlim([-3.5432    1.8662]);
ylim([ -5.6786    2.5062]);
print(gcf,'netscience-figure-cut.eps','-depsc2','-painters')
%% Compute 1-norm
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
z = x(2:end-1); 
example_netscience_fig(A,xy,max(ns*z,nm),sqrt(z))
title(sprintf('%i nonzeros',nnz(z > 1e-5)));
xlim([-3.5432    1.8662]);
ylim([ -5.6786    2.5062]);
print(gcf,'netscience-figure-acl.eps','-depsc2','-painters')

%%
!cp netscience-figure-* ~/Dropbox/publications/alg-anti-diff-icml/anti-diff-paper-1/