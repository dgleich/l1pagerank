function [x,result] = prl1_gurobi(A,alpha,s,d,kappa)
% PRCUT_GUROBI Solve the PR cut problem using Gurobi

n = size(A,1);

Q = alpha*diag(d) + (diag(d) - A);
%cvx_begin
%    variable x(n);
%    minimize 1/2*quad_form(x,Q) - alpha*sum(x.*d.*s) + kappa*sum(d.*x)
%    subject to
%        x >= 0;
%cvx_end
%yacl = x;

clear model;
model.obj = full(kappa*d - alpha*d.*s);
model.Q = sparse(0.5*Q);
model.lb = zeros(n,1);
model.ub = Inf*ones(n,1);
model.A = sparse(0,n);
model.rhs = [];
model.sense = '';

clear params;
params.outputflag = 0;
params.OptimalityTol = sqrt(eps(1));

result = gurobi(model,params);
x = result.x; % extract the solution


            