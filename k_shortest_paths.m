classdef k_shortest_paths < handle
    %K_SHORTEST_PATHS Finds the k shortest paths in a graph.

    properties (GetAccess = public, SetAccess = protected)
        %G The graph
        G;
        
        %V_SOURCE The source vertex
        v_source;
        
        %V_SINK The sink vertex
        v_sink;
                
        %LAST_PATHS A cell array containing all shortest paths found so far
        last_paths;
        
        %LAST_COSTS A vector storing the costs of all shortest paths found so far
        last_costs;
        
        %SHORTEST_PATH_FUN The function to use for finding the shortest paths
        shortest_path_fun;
        
        %FIRST_SHORTEST_PATH_FUN The function to use for finding the first shortest path
        first_shortest_path_fun;
    end
    
    properties (GetAccess = public, SetAccess = private)        
        %ITERATION The last iteration
        iteration = 0;
        
        %PATH_NOT_FOUND True if no path was found in the last iteration
        path_not_found = false;
        
        %VISITOR The visitor object, if any
        visitor;
    end
    
    methods 
        function set_visitor(obj, visitor)
        %SET_VISITOR Sets the visitor object.
        %   set_visitor(visitor) sets visitor as the new visitor of this object.
        %   visitor must be of type k_shortest_paths_visitor or empty, if you want
        %   to clear the visitor.
            
            assert(isempty(visitor) || isa(visitor, 'k_shortest_paths_visitor'), 'visitor must be of type k_shortest_paths_visitor.');
            
            obj.visitor = visitor;
        end
        
        function [ paths, costs ] = find(obj, k)
        %FIND Finds shortest paths in the graph.
        %   [paths, costs] = obj.find(k) finds the k shortest paths in the
        %   graph. Alternatively, you can use this function to iteratively find
        %   shortest paths, based on some condition. In this case, call find()
        %   without specifying k. The first call to find() will compute the
        %   shortest path, the second call the two shortest paths and so on.
            
            if nargin >= 2
                assert(k > 0, 'k must be greater 0.');
                assert(obj.iteration == 0, 'k is only allowed at the first call to find().');
            end
            
            if nargin >= 2
                for i = 1:k
                    [ paths, costs ] = obj.find();
                end
            else
                [ paths, costs ] = obj.next_iteration();
                
                if ~isempty(obj.visitor)
                    [ paths, costs ] = obj.visitor.paths_computed(obj, paths, costs);
                end
                
                obj.iteration = obj.iteration + 1;
            end
        end
    end
    
    methods (Access = protected)
        function on_path_not_found(obj)
            obj.path_not_found = true;
            warning('k_shortest_paths:path_not_found', 'No shortest path found.');
        end
    end
    
    methods (Abstract = true, Access = protected)
        [ paths, costs ] = next_iteration(obj);
    end
    
    methods (Static = true, Access = protected)
        function input_parser = get_ctor_input_parser()
            persistent parser;
            
            if isempty(parser)
                parser = inputParser();
                parser.addRequired('G', @(x) ~isempty(x) && size(x,2) == 3);
                parser.addRequired('v_source', @isscalar);
                parser.addRequired('v_sink', @isscalar);
                parser.addRequired('shortest_path_fun', @(x) isa(x, 'function_handle'));
                parser.addOptional('first_shortest_path_fun', [], @(x) isa(x, 'function_handle'));
            end
            
            input_parser = parser;
        end
        
        function arcs = find_common_arcs(path_a, path_b) 
        %FIND_COMMON_ARCS Finds all arcs that are part of two paths.
        %   arcs = find_common_arcs(path_a, path_b) finds all arcs that are path
        %   of both path_a and path_b. It is assumed It is assumed that both parts
        %   have the same source and sink, and that there are no cycles. arcs is a
        %   n*2 matrix describing all common arcs or empty if no such arcs were found.
            
            arcs = [];
        
            common_vertices = intersect(path_a, path_b);
            common_vertices = unique(common_vertices);
            common_vertices = setdiff(common_vertices, [path_a(1) path_a(end)]);
            
            if isempty(common_vertices)
                % no common vertices, thus no common arcs
                return
            end
            
            for i = 1:numel(common_vertices)
                from = common_vertices(i);
                
                idx_path_a = find(path_a == from, 1);
                idx_path_b = find(path_b == from, 1);
                
                to_path_a = path_a(idx_path_a + 1);
                to_path_b = path_b(idx_path_b + 1);
                
                if to_path_a == to_path_b
                    arcs = vertcat(arcs, [from to_path_a]);
                end
            end
        end
        
        function graph = paths_to_graph(G, paths)
        %PATHS_TO_GRAPH Creates a graph from a number of paths.
        %   graph = paths_to_graph(G, paths) creates a graph graph which is a
        %   n*3 matrix of vertices and their weights from a number of paths in
        %   G, which is also a n*3 matrix. paths is a vector describing a path or
        %   a cell array of paths.
            
            arc_idx = 1;
            G_arcs = [G(:,1) G(:,2)];
            
            if iscell(paths)
                num_arcs = sum(cellfun(@(x) numel(x) - 1, paths));
                graph = zeros(num_arcs, 3);
                cellfun(@path_to_arcs, paths);
            else
                num_arcs = numel(paths) - 1;
                graph = zeros(num_arcs, 3);
                path_to_arcs(paths);
            end

            function path_to_arcs(path)
                for i = 1:numel(path) - 1
                    from = path(i);
                    to = path(i + 1);
                    
                    [ ~, idx ] = ismember([from to], G_arcs, 'rows');
                    
                    graph(arc_idx,:) = [from to G(idx,3)];
                    arc_idx = arc_idx + 1;
                end
            end
        end
        
        function graph = generate_combined_graph(graph_1, graph_2)
        %GENERATE_COMBINED_GRAPH Generates a combined graph from two graphs
        %   graph = generate_combined_graph(graph_1, graph_2) generates a
        %   combined graph graph which is a n*3 matrix of all arcs and their
        %   weights. graph_1 and graph_2 are assumed to be of identical
        %   structure. graph contains all vertices and arcs of graph_1 and
        %   graph_2, but with multi-arcs and arcs that go in both directions (a
        %   <-> b) removed.
            
            graph = vertcat(graph_1, graph_2);
            graph = unique(graph, 'rows');
            
            arcs = graph(:,1:2);

            % we need to check for arcs a <-> b in the combined graph, i.e.,
            % arcs in both directions, and remove them. to do this we conatenate
            % arcs with the flipped version of arcs. If there are no connections
            % in both directions, there are no duplicate rows in the resulting
            % matrix. All non-unique rows represent arcs in both directions.
            
            arcs = vertcat(arcs, fliplr(arcs));            
            [~, unique_idx, ~] = unique(arcs, 'rows', 'last');
            bidir_idx = setdiff(1:size(arcs, 1), unique_idx);
            
            graph(bidir_idx,:) = [];
        end
        
        function [ paths, costs ] = find_paths_in_combined_graph(graph, v_source, v_sink)
        %FIND_PATHS_IN_COMBINED_GRAPH Finds all paths in a combined graph.
        %   [ paths, costs ] = find_paths_in_combined_graph(graph, v_source, v_sink)
        %   finds all paths from v_source to v_sink in the combined graph graph,
        %   which is a n*3 vector of all arcs and their weights. a combined
        %   graph is a graph created with the function generate_combined_graph,
        %   i.e., a DAG. paths is a cell array of all paths found, costs is a
        %   vector that describes the total weights of the paths.
        
            v_max = max(max(graph(:,1:2)));
            
            % get all arcs from v_source
        
            source_arcs_idx = find(graph(:,1) == v_source);
            num_source_arcs = size(source_arcs_idx, 1);
            
            paths = cell(1, num_source_arcs);
            costs = zeros(1, num_source_arcs);
            
            % for every arc from v_source, compute the shortest path to v_sink
            % (which is fast, since the graph is acyclic). Then, remove all arcs
            % that comprise the shortest path and start again. This works since
            % all shortest paths are arc disjoint.
            
            for i = 1:num_source_arcs
                [ p, c ] = dag_shortest_paths(graph, v_source);
                
                path = p{v_sink};
                cost = c(v_sink);


                if isempty(path)
                    error('v_sink unreachable from v_source, invalid graph.');
                end
                
                paths{i} = path;
                costs(i) = cost;
                
                arcs = k_shortest_paths.path_to_arcs(path);
                [ ~ , arcs_idx ] = ismember(arcs, graph(:,1:2), 'rows');
                
                graph(arcs_idx, :) = [];
            end
        end
        
        function R = reverse_arcs(G, paths)
        %REVERSE_ARCS Reverses all arcs in G that are part of any path in paths.
        %   R reverses all arcs and negates their weight of
        %   all arcs that are part of any path in paths. paths is either a
        %   vector [v1 v2 ... vn] describing one path in G or a cell array of
        %   paths. All arcs are assumed to be unique among all paths. G is a
        %   n*3 matrix. R is the modified graph, also n*3.
            
            G_arcs = [G(:,1) G(:,2)];
            
            R = G;
            R_remove = [];
            R_add = [];
            
            if iscell(paths)
                cellfun(@reverse_arcs, paths);
            else
                reverse_arcs(paths);
            end
            
            R_add = unique(R_add, 'rows');
            
            R = setdiff(R, R_remove, 'rows');
            R = vertcat(R, R_add);
       
            function reverse_arcs(path)
                for v = 1:numel(path) - 1
                    from = path(v);
                    to = path(v + 1);
                    
                    [ ~, idx ] = ismember([from to], G_arcs, 'rows');
                    
                    weight = G(idx,3);
                    
                    R_remove = vertcat(R_remove, [from to weight]);
                    R_add = vertcat(R_add, [to from -weight]);
                end  
            end
        end
        
        function arcs = path_to_arcs(path)
        %PATH_TO_ARCS Returns a n*2 matrix of all arcs in path.
            
            arcs = zeros(numel(path) - 1, 2);
            
            for i = 1:numel(path) - 1
                arcs(i,:) = [path(i) path(i + 1)];
            end
        end
        
        function weight = arc_weight(G, arc)
        %ARC_WEIGHT Returns the weight of arc [a b] in G, where G is a n*3 
        %   matrix [from to weight ; ...].
            
            weight = inf;
            
            G_arcs = [G(:,1) G(:,2)];
            [ ~, idx ] = ismember([arc(1) arc(2)], G_arcs, 'rows');
            
            if idx > 0
                weight = G(idx,3);
            end
        end
    end
end
