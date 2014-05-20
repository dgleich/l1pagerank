Algorithm Anti-Differentiation: A case study with PageRank, min-cuts, and flow
================

### David F. Gleich
### Michael W. Mahoney

Contents
--------

These are the computational codes and data that go along with the ICML 2014 paper

    David F. Gleich and Michael Mahoney
    Algorithm Anti-Differentiation: A case study with min-cuts, spectral, and flow
    In Proceedings of the International Conference on Machine Learning, 2014.

Files
-----

### Demos

* `figure_1_example.m` Produce the vectors for Figure 1 in the ICML paper
* `example_netscience.m` Produce the netscience Figures for the ICML paper
* `demos/simple_example.m` A simple example suitable for live demos

### Algorithms

* `acl_method.m` An implementation of the Andersen-Chung-Lang
          push method to compute a PageRank vector
          
* `flow_improve.m` An algorithm for the FlowImprove method of Andersen and Lang

* `prl1_gurobi.m` Use gurobi to solve an l1-PageRank problem
* `prcut_gurobi.m` Use gurobi to solve a min-cut on the PageRank cut graph

### Helpers

* `cut_graph.m` Construct the cut-graph for FlowImiprove
* `general_cut_graph.m` Construct the PageRank cut-graph for general constructions.
* `gmatrices.m` Return a struct of all the matrices associated with a graph.
* `igraph_draw.m` Use igraph via python to layout a graph
* `incidence_matrix.m` Create the node-edge indicence matrix for a graph
* `laplacian.m` Consruct the Laplacian matrix of a graph
* `load_graph.m` Load one of the datasets and avoid icky path problems
* `nlaplacian.m` Construct the normalized Laplacian of a graph
* `normout.m` Compute a row-stochastic matrix from an adjacency matrix.
* `set_figure_size.m` The most useful script to make Matlab figures look good
* `setup_paths.m` Add all the required libraries to the matlab path, and the code itself
* `wadjacency.m` Construct the weighted adjacency matrix of a graph

### Testing

* `test_prcut_gurobi.m` Make sure gurobi and CVX agree for the prcut problem
* `test_prl1_gurobi.m` Make sure gurobi and CVX agree for the prl1 problem
          
          
### Exploration

This directory contains a whole bunch of scripts that we developed while
working out the theory. Think of this as a lab-notebook for what we looked
at before the full picture came together.
