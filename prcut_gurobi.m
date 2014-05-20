function [x,result] = prcut_gurobi(A,alpha,s,d)
% PRCUT_GUROBI Solve the PR cut problem using Gurobi

sbar = d - s;
AS = [0              alpha*(s)'   0;
      alpha*(s)  A     alpha*(sbar);
      0              alpha*(sbar)' 0];
n = size(AS,1);
      
[ei ej ev] = find(triu(AS,0));
C = diag(sparse(ev));
nedge = numel(ei);
Bu = sparse([(1:nedge)'; (1:nedge)'], [ei ej], ...
            [ones(nedge,1); -1*ones(nedge,1)], ...
            nedge, n);
            
m = nedge;
c = zeros(m + n, 1);
c(1:nedge) = 1;
lb = zeros(m + n, 1);
ub = Inf*ones(m + n, 1);

% Constraint matrix
% Matrix*x <= rhs
e1 = zeros(m+n,1)'; e1(m+1) = 1;
en = zeros(m+n,1)'; en(m+n) = 1;
Matrix = [e1 % x(1) = 1
          en % x(n) = 0
          -speye(m) -C*Bu % -y <= C*Bu*x 
          -speye(m) C*Bu  % C*Bu*x <= y ----> -y + C*Bu*x <= 0 
          ];
rhs = [1; 0; zeros(2*m,1)];

clear model; 
model.obj = c;
model.A = sparse(Matrix);
model.rhs = rhs;
model.sense = ['=';'=';repmat('<',2*m,1)];
model.lb = lb;
model.ub = ub;
model.vtype = 'C';

clear params;
params.outputflag = 0;
%params.OptimalityTol = sqrt(eps(1));

result = gurobi(model,params);
x = result.x(m+2:end-1); % extract the solution


            