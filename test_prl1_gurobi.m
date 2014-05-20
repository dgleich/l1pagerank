figure_1_example; % run things with cvx
%%
yg = prl1_gurobi(A,alpha,Sv,M.d,kappa);
assert(all(abs(yr - yg) < 10*sqrt(eps(1))));
