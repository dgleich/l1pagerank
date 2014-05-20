function h = example_netscience_fig(A,xy,s,c)

h = [];
[px,py] = gplot(A,xy);
h(end+1) = plot(px,py,'k-','LineWidth',0.25); hold on;
f = c.^2 > 1e-5;
h(end+1) = scatter(xy(:,1),xy(:,2),s,c,'Filled');
h(end+1) = scatter(xy(~f,1),xy(~f,2),s(~f));
set(h(end),'MarkerEdgeColor','k');
colormap(flipud(hot))
hold off;
axis off;
set_figure_size([2.5,2.5]);