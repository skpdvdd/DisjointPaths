function [ paths, costs ] = dag_shortest_paths( G, v_source )
%DAG_SHORTEST_PATHS Finds the shortest paths from v_source in a graph.
%   [paths, costs] = dag_shortest_paths( G, v_source ) finds the
%   shortest path from v_source to all other vertices in a DAG G (must be a n*3
%   matrix [from to weight ; ...] with positive weights) and returns the paths
%   and costs. paths{i} contains the shortest path from v_source to vertex i.
%   costs contains the path costs (Inf for unreachable vertices).

    G = sortrows(G, [2 1]);
    v = G(:,3);
    G(:,3) = 1;
    v_max = max([G(:,1) ; G(:,2)]);
    
    G = sparse(G(:,1), G(:,2), G(:,3), v_max, v_max);

    [ costs, paths, ~ ] = graphshortestpath(G, v_source, 'Method', 'Acyclic', 'Weights', v);
end
