classdef k_shortest_vertex_disjoint_paths < k_shortest_paths
    %K_SHORTEST_VERTEX_DISJOINT_PATHS Finds the k shortest vertex disjoint paths in a graph.

    methods
        function obj = k_shortest_vertex_disjoint_paths(G, v_source, v_sink, shortest_path_fun, varargin)
        %K_SHORTEST_VERTEX_DISJOINT_PATHS Creates a new object of this class.
        %   k_shortest_vertex_disjoint_paths(G, v_source, v_sink, shortest_path_fun, first_shortest_path_fun)
        %   creates a new object for finding vertex disjoint shortest paths in G (a n*3
        %   matrix [from to weight ; ...]) from v_source to v_sink. shortest_path_fun is a handle to
        %   the function to use for computing shortest paths.
        %   first_shortest_path_fun (optional) is a handle to the function to
        %   use when computing the first shortest path.
            
            parser = k_shortest_paths.get_ctor_input_parser();
            parser.parse(G, v_source, v_sink, shortest_path_fun, varargin{:});            
            
            obj.G = G;
            obj.v_source = v_source;
            obj.v_sink = v_sink;
            obj.shortest_path_fun = fcnchk(shortest_path_fun);
            
            if ~isempty(varargin)
                obj.first_shortest_path_fun = fcnchk(varargin{1});
            end
        end
    end
    
    methods (Static = true, Access = protected)
        function [ result, v_offset ] = add_concomitant_vertices(graph, paths, original_graph)
            
            % add concomitant vertices: Let every path in paths be v1 ... vk. For
            % every vi with i=2:k-1 add a concomitant vertex vi' and create
            % an arc to vi-1, with the weight being the negated weight of
            % the arc from vi-1 to vi
            
            % the concomitant vertices need to have ids that can be associated with
            % the original vertices, i.e. vi <-> vi'. Thus, we compute the
            % max vertex id and set the id of vi' to be vi + v_max
            
            v_max = max(max(graph(:,1:2)));
            v_source = paths{1}(1);
            v_sink = paths{1}(end);
            v_offset = v_max;
            
            path_arcs = [];
            concomitant_vertices = [];
            
            if iscell(paths)
                cellfun(@add_concomitant_arcs, paths);
            else
                add_concomitant_arcs(paths);
            end
            
            function add_concomitant_arcs(path)
                if numel(path) == 2
                    return
                end
                
                arcs = k_shortest_paths.path_to_arcs(path);
                
                path_arcs = vertcat(path_arcs, arcs);
                
                arcs = arcs(2:end-1,:);
                num_arcs = size(arcs, 1);
                
                for j = 1:num_arcs
                    from = arcs(j,1);
                    to = arcs(j,2);
                    weight = k_shortest_paths.arc_weight(original_graph, [from to]);
                    
                    concomitant_vertex = [to + v_max , from , -weight];
%                     fprintf('new concomit vertex %s\n', mat2str(concomitant_vertex));
                    
                    concomitant_vertices = vertcat(concomitant_vertices, concomitant_vertex);
                end
            end
            
            path_arcs = unique(path_arcs, 'rows');
            concomitant_vertices = unique(concomitant_vertices, 'rows');
            
            graph = vertcat(graph, concomitant_vertices);
            
            % find all arcs in the original graph that are not part of any
            % shortest path. check these arcs for arcs to vertices b that are
            % part of the shortest paths and redirect to b'.
                        
            original_graph_arcs = original_graph(:,1:2);
            original_graph_arcs_no_path = setdiff(original_graph_arcs, path_arcs, 'rows');
            
            % using a lookup table for fast check if b is in any a->b arc of any
            % shortest path. we do not replace vertices pointing to the source
            % or sink, so we set the lut to false explicitly.
            
            vertex_lut = false(1, max(max(original_graph_arcs)));
            
            for i = 1:size(path_arcs, 1)
                vertex_lut(path_arcs(i,2)) = true;
            end
            
            vertex_lut(v_source) = false;
            vertex_lut(v_sink) = false;
            
            graph_arcs = graph(:,1:2);
            
            for i = 1:size(original_graph_arcs_no_path)
                from = original_graph_arcs_no_path(i,1);
                to = original_graph_arcs_no_path(i,2);
                
                if vertex_lut(to) == true
                    % if there is an arc from a -> b where b is a in vertex
                    % in any shortest path, change to a -> b'

                    [ ~, graph_idx ] = ismember([from to], graph_arcs, 'rows');

                    graph(graph_idx,2) = to + v_max;
                end
            end
            
            result = graph;
        end
    end
        
    methods (Access = protected)
        function [ paths, costs ] = next_iteration(obj)
        %NEXT_ITERATION Computes the next iteration, i.e. the i+1 shortest paths.
        %   [paths, costs] computes the k+1 shortest paths (where k was the
        %   iteration before calling this method) and returns them. paths is a
        %   cell array of size i+1, with every element being a vector that
        %   describes a shortest path [v1 v2 ... vn]. costs is a vector of size
        %   i+1 holding the cost for every path.
            
            G = obj.G;
        
            if ~isempty(obj.visitor)
                G = obj.visitor.begin(obj, G);
            end
        
            if obj.iteration == 0
                fun = obj.shortest_path_fun;

                if ~isempty(obj.first_shortest_path_fun)
                    fun = obj.first_shortest_path_fun;
                end

                [ p, c ] = fun(G, obj.v_source);
                
                if ~isempty(obj.visitor)
                    [ p, c ] = obj.visitor.shortest_paths_computed(obj, p, c);
                end
                
                path = p{obj.v_sink};
                costs = c(obj.v_sink);
                
                if isempty(path)
                    obj.on_path_not_found();
                    paths = [];
                    costs = [];
                    return
                end
                
                obj.last_paths = cell(1);
                obj.last_paths{1} = path;
                obj.last_costs = costs;
                
                paths = obj.last_paths;
            else
                if obj.path_not_found
                    error('k_shortest_paths:last_shortest_path_not_found', 'No shortest path found at last iteration.');
                end
                
                % reverse all arcs that are in any shortest path to obtain R
                        
                R = k_shortest_paths.reverse_arcs(G, obj.last_paths); 
                
                if ~isempty(obj.visitor)
                    R = obj.visitor.graph_reversed(obj, R);
                end
                
                % add concomitant vertices and redirect arcs
                
                [ R, v_offset ] = k_shortest_vertex_disjoint_paths.add_concomitant_vertices(R, obj.last_paths, G);
                
                if ~isempty(obj.visitor)
                    [ R , v_offset ] = obj.visitor.concomitant_vertices_added(obj, R, v_offset);
                end
                
                % find the shortest path in the transformed graph and find all
                % vertices common to the found path and the shortest paths
                
                [ p, c ] = obj.shortest_path_fun(R, obj.v_source);
                
                if ~isempty(obj.visitor)
                    [ p, c ] = obj.visitor.shortest_paths_computed(obj, p, c);
                end
                
                path = p{obj.v_sink};
                cost = c(obj.v_sink);
                
                if isempty(path)
                    obj.on_path_not_found();
                    paths = [];
                    costs = [];
                    return
                end
                
                paths_vertices = horzcat(obj.last_paths{:});
                paths_vertices = unique(paths_vertices);
                
                common_vertices = intersect(path, paths_vertices);
                common_vertices = setdiff(common_vertices, [obj.v_source obj.v_sink]);
                
                if isempty(common_vertices)
                    % if there are no common vertices, we are done
                    
                    obj.last_paths{end+1} = path;
                    obj.last_costs(end+1) = cost;
                    
                    paths = obj.last_paths;
                    costs = obj.last_costs;
                    
                    return
                end
                
                %----

                % otherwise combine all vertices v in path with their corresponding
                % concomitant vertices v' ... that is change v' to v. the relation
                % between the ids is v' = v + v_offset
                
                % find all concomitant vertices in the found path and transform
                % them to their "original" vertices
                
                arcs_combined = k_shortest_paths.path_to_arcs(path);
                
                for i = 1:size(arcs_combined, 1)
                    from = arcs_combined(i,1);
                    to = arcs_combined(i,2);
                    
                    arcs_combined(i,3) = k_shortest_paths.arc_weight(R, [from to]);
                    
                    if from > v_offset, arcs_combined(i,1) = from - v_offset; end
                    if to > v_offset, arcs_combined(i,2) = to - v_offset; end
                end
                                    
                % generate the combined graph
                
                shortest_path_arcs = k_shortest_paths.paths_to_graph(G, obj.last_paths);
                
                combined_graph = k_shortest_paths.generate_combined_graph(shortest_path_arcs, arcs_combined);
                
                if ~isempty(obj.visitor)
                    combined_graph = obj.visitor.combined_graph_generated(obj, combined_graph);
                end
                                
                % the combined graph holds all shortest path

                [obj.last_paths, obj.last_costs] = k_shortest_paths.find_paths_in_combined_graph(combined_graph, obj.v_source, obj.v_sink);

                paths = obj.last_paths;
                costs = obj.last_costs;
            end
        end
    end 
end

