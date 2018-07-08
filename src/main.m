function main(png_filename)
V = [];
E = [];

if (nargin < 1)
%     nX = 3;
%     nY= 3;
%     
%     [X,Y] = meshgrid(linspace(-5,5,nX), linspace(-2.5,2.5,nY));
%     V = [X(:), Y(:)];
%     
%     %build trusses by triangulating nodes
%     F = delaunayTriangulation(V(:,1), V(:,2));
%     E = edges(F);
    
    %IMMA DRAW THE T
    V = [0 2;
         0 3;
         1 3;
         2 3;
         3 3;
         3 2;
         2 2;
         2 1;
         2 0;
         1 0;
         1 1;
         1 2];
          
    E = [1 2;
         2 3;
         3 4;
         4 5;
         5 6;
         6 7;
         7 8;
         8 9;
         9 10;
         10 11;
         11 12;
         12 1];
    
else
    [V, E] = bwmesh(png_filename);
end


%TODO: find a better way to do this. BFS fails if the edges are not unique
%uniquify
%E = E(E(:, 1) < E(:, 2), :);

figure
hold on
line([V(E(:,1),1)';V(E(:,2),1)'],[V(E(:,1),2)';V(E(:,2),2)'], 'Color', [0 0 1]);
hold off

levelManager = leveller(V, E);

for i = 1:levelManager.maxLevel
    Einds = levelManager.getEdgesBelowLevel(i);
    levelE = E(Einds, :);
    levelV = levelManager.getVerticesBelowLevel(i);
    
    tensions = topologyOptimization(V, levelE, levelV, i);
end

end