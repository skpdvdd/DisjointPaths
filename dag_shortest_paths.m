function [ paths, costs ] = dag_shortest_paths( G, v_source )
%DAG_SHORTEST_PATHS Finds the shortest paths from v_source in a graph.
%   [paths, costs] = dag_shortest_paths( G, v_source ) finds the
%   shortest path from v_source to all other vertices in a DAG G (must be a n*n
%   sparse matrix) and returns the paths and costs. paths{i} contains the shortest
%   path from v_source to vertex i. costs contains the path costs (Inf for
%   unreachable vertices).

    [ costs, paths, ~ ] = graphshortestpath(G, v_source, 'Method', 'Acyclic');
end
