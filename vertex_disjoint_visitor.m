classdef vertex_disjoint_visitor < handle
    %VERTEX_DISJOINT_VISITOR Visitor class for use with k_shortest_vertex_disjoint_paths.
    
    methods
        R = graph_reversed(obj, subject, R);
        %GRAPH_REVERSED Called after the graph was reversed.
        %   R is the resulting graph, a sparse matrix.
        
        graph = combined_graph_generated(obj, subject, graph);
        %COMBINED_GRAPH_GENERATED Called after the combined graph was generated.
        %   graph is the combined graph, a n*3 matrix [from to weight ; ...].
        
        [ graph , v_offset ] = concomitant_vertices_added(obj, subject, graph, v_offset);
        %CONCOMITANT_VERTICES_ADDED Called after the reversed graph was
        %   augmented with the concomitant vertices. graph is the resulting
        %   graph, a n*3 matrix, v_offset is the offset between the vertices and
        %   the corresponding concomitant vertices: v' = v + v_offset.
    end 
end
