classdef k_shortest_arc_disjoint_paths < k_shortest_paths
    %K_SHORTEST_ARC_DISJOINT_PATHS Finds the k shortest arc disjoint paths in a graph.
        
    methods
        function obj = k_shortest_arc_disjoint_paths(G, v_source, v_sink, shortest_path_fun, varargin)
        %K_SHORTEST_ARC_DISJOINT_PATHS Creates a new object of this class.
        %   k_shortest_arc_disjoint_paths(G, v_source, v_sink, shortest_path_fun, first_shortest_path_fun)
        %   creates a new object for finding arc disjoint shortest paths in G (a n*3
        %   matrix [from to weight ; ...] with positive weights) from v_source to v_sink.
        %   shortest_path_fun is a handle to the function to use for computing shortest paths.
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
    
    methods (Access = protected)
        function [ paths, costs ] = next_iteration(obj)  
        %NEXT_ITERATION Computes the next iteration, i.e. the i+1 shortest paths.
        %   [paths, costs] = next_iteration() computes the k+1 shortest paths (where k was the
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
                
                R = k_shortest_paths.reverse_arcs(G, obj.last_paths); 
                
                if ~isempty(obj.visitor)
                    R = obj.visitor.graph_reversed(obj, R);
                end
                
                [ p, c ] = obj.shortest_path_fun(R, obj.v_source);
                
                if ~isempty(obj.visitor)
                    [ p, c ] = obj.visitor.shortest_paths_computed(obj, p, c);
                end
                
                path = p{obj.v_sink};
                cost = c(obj.v_sink);
                
                if isempty(path)
                    obj.on_path_not_found();
                    paths = [];
                    costs = Inf;
                    return
                end
                
                common_arcs = k_shortest_paths.find_common_arcs(obj.last_paths{1}, fliplr(path));
                
                if isempty(common_arcs)
                    obj.last_paths{end + 1} = path;
                    obj.last_costs(end + 1) = cost;
                    
                    paths = obj.last_paths;
                    costs = obj.last_costs;
                else
                    arcs_1 = k_shortest_paths.paths_to_graph(G, obj.last_paths);
                    arcs_2 = k_shortest_paths.paths_to_graph(R, path);
                                        
                    graph = k_shortest_paths.generate_combined_graph(arcs_1, arcs_2);
                    
                    if ~isempty(obj.visitor)
                        graph = obj.visitor.combined_graph_generated(obj, graph);
                    end
                              
                    [obj.last_paths, obj.last_costs] = k_shortest_paths.find_paths_in_combined_graph(graph, obj.v_source, obj.v_sink);

                    paths = obj.last_paths;
                    costs = obj.last_costs;
                end
            end
        end
    end 
end
