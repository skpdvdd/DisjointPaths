function [ path, cost ] = bellman_ford_shortest_paths( G, v_source, v_sink )
%BELLMAN_FORD_SHORTEST_PATHS Finds the shortest path in a graph.
%   [path, cost] = bellman_ford_shortest_paths( G, v_source, v_sink ) finds the
%   shortest path from v_source to v_sink in a graph G (must be a sparse matrix)
%   and returns the path and the total cost. If no path could be found, path is
%   empty and cost is Inf.

    [ cost, path, ~ ] = graphshortestpath(G, v_source, v_sink, 'Method', 'Bellman-Ford');
end

