classdef arc_disjoint_visitor < handle
    %ARC_DISJOINT_VISITOR Visitor class for use with k_shortest_arc_disjoint_paths.
    
    methods (Abstract = true)
        R = graph_reversed(obj, subject, R);
        %GRAPH_REVERSED Called after the graph was reversed.
        %   R is the resulting graph, a sparse matrix.
        
        graph = combined_graph_generated(obj, subject, graph);
        %COMBINED_GRAPH_GENERATED Called after the combined graph was generated.
        %   graph is the combined graph, a n*3 matrix [from to weight ; ...].
    end
end
