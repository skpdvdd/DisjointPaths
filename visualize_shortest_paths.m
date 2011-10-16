function visualize_shortest_paths(G, paths)
%VISUALIZE_SHORTEST_PATHS Creates a visualization a graph and its shortest paths.
%   visualize_shortest_paths(G, paths) creates a visualization of the graph G
%   and its shortest paths. G must be a quadratic sparse matrix and paths must
%   be a cell array with k elements, with each one being a vector describing
%   a path in the form [v1 v2 ... vn].

    v_max = max([G(:,1) ; G(:,2)]);
    v_0 = find(G(:,3) == 0);
    G(G(:,3)==0,3) = -Inf;
    G = sparse(G(:,1), G(:,2), G(:,3), v_max, v_max);

    h = view(biograph(G, [], 'ShowWeights', 'on'));
    
    hue = 0;
    hue_step = 0;
    
    if numel(paths) > 1
        hue_step = 1 / (numel(paths) - 1);
    end
    
    for i = 1:numel(paths)
        color = hsv2rgb([hue 0.9 0.9]);
        
        edges = getedgesbynodeid(h, get(h.Nodes(paths{i}), 'ID'));
        set(edges, 'LineColor', color);
        set(edges, 'LineWidth', 1.5);
        
        hue = hue + hue_step;
    end
end
