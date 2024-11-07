%% 
function [pos,n] = find_closest(G,x)
    [n,pos] = min(vecnorm(G.cells.centroids -x,2,2));
end