function A = general_cut_graph(G,ws,wsbar)
ws = ws(:);
wsbar = wsbar(:);
A = [0 ws' 0; ws G wsbar; 0 wsbar' 0];
 
