classdef arc_disjoint_visitor_example < arc_disjoint_visitor
    %ARC_DISJOINT_VISITOR_EXAMPLE Visitor example.

    methods
        function R = graph_reversed(obj, subject, R)
            fprintf('reversed graph of iteration %d:\n', subject.iteration);
            full(R)
        end
        
        function graph = combined_graph_generated(obj, subject, graph)
            fprintf('combined graph of iteration %d:\n', subject.iteration);
            graph
        end
    end   
end
