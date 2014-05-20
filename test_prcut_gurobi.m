figure_1_example; % run things with cvx
%%
[xc,res] = prcut_gurobi(G,alpha,M.d.*Sv,M.d);
assert(all(abs(zc - xc) < 10*eps(1)));