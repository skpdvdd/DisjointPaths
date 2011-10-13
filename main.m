
clear all;
clc;

child_handles = allchild(0);
names = get(child_handles,'Name');
k = find(strncmp('Biograph Viewer', names, 15));
close(child_handles(k));

% -----------------------------

G = [1 2 0.001
     1 3 0.001
     1 4 0.001
     2 5 -2.8
     2 6 -2.8
     3 5 2.5
     3 6 2.5
     3 7 2.5
     4 6 -2.7
     4 7 -2.7
     5 8 2.9
     5 9 2.9
     6 8 -2.5
     6 9 -2.5
     6 10 -2.5
     7 9 -2.6
     7 10 -2.6
     8 11 -2.9
     9 11 -2.85
     10 11 2.2];
     
G = sparse(G(:,1), G(:,2), G(:,3), 11, 11);

% kp = k_shortest_arc_disjoint_paths(G, 1, 11, @bellman_ford_shortest_paths);

kp = k_shortest_vertex_disjoint_paths(G, 1, 11, @bellman_ford_shortest_paths, @dag_shortest_paths);
 
% -------------------------------

% G = [1 2 2
%      2 3 2
%      3 4 2
%      4 6 2
%      1 4 7
%      2 5 4
%      5 6 4];
%  
% G = sparse(G(:,1), G(:,2), G(:,3), 6, 6);
% 
% kp = k_shortest_arc_disjoint_paths(G, 1, 6, @bellman_ford_shortest_paths);
% kp = k_shortest_vertex_disjoint_paths(G, 1, 6, @bellman_ford_shortest_paths);

% --------------------------------

% G = [1 2 1
%      2 3 1
%      3 4 1
%      4 5 1
%      5 6 1
%      6 7 1
%      1 4 5
%      4 6 3
%      2 5 4
%      5 7 3];
%  
% G = sparse(G(:,1), G(:,2), G(:,3), 7, 7);
% 
% kp = k_shortest_arc_disjoint_paths(G, 1, 7, @bellman_ford_shortest_paths);
% kp.set_visitor(arc_disjoint_visitor_example());
% kp = k_shortest_vertex_disjoint_paths(G, 1, 7, @bellman_ford_shortest_paths);

% --------------------------------

k = 5;
vis_all = true;

for i = 1:k
    [ paths, costs ] = kp.find();
    
    if isempty(paths)
        disp('No more shortest paths found.');
        break
    end
    
    fprintf('Found %d shortest paths:\n', i);
    for j = 1:numel(paths)
        fprintf('  %s (cost: %.3f)\n', mat2str(paths{j}), costs(j));
    end
    
    if vis_all
        visualize_shortest_paths(G, paths);
        input('Press any key to continue..');
    end
    
    disp('------');
end
