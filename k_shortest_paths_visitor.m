classdef k_shortest_paths_visitor < handle
    %K_SHORTEST_PATHS_VISITOR Visitor class for use with k_shortest_vertex_*_paths.
    
    methods
        R = graph_reversed(obj, subject, R);
        %GRAPH_REVERSED Called after the graph was reversed.
        %   R is the resulting graph, a sparse matrix.
        
        graph = combined_graph_generated(obj, subject, graph);
        %COMBINED_GRAPH_GENERATED Called after the combined graph was generated.
        %   graph is the combined graph, a n*3 matrix [from to weight ; ...].
        
        [ paths, costs ] = shortest_paths_computed(obj, subject, paths, costs);
        %SHORTEST_PATHS_COMPUTED Called after computing the shortest paths in G (or a
        %   transformed version of G) from v_source to all other vertices. paths
        %   is a cell array of shortest paths from v_source to all other
        %   vertices, costs a vector of path costs.
        
        [ graph , v_offset ] = concomitant_vertices_added(obj, subject, graph, v_offset);
        %CONCOMITANT_VERTICES_ADDED Called after the reversed graph was
        %   augmented with the concomitant vertices. graph is the resulting
        %   graph, a n*3 matrix, v_offset is the offset between the vertices and
        %   the corresponding concomitant vertices: v' = v + v_offset.
    end 
end
