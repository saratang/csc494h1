function main(png_filename)
V = [];
E = [];

if (nargin < 1)
    nX = 3;
    nY= 4;
    
    [X,Y] = meshgrid(linspace(-5,5,nX), linspace(-2.5,2.5,nY));
    V = [X(:), Y(:)];
    
    %build trusses by triangulating nodes
    F = delaunayTriangulation(V(:,1), V(:,2));
    E = edges(F);
    
else
    [V, E] = bwmesh(png_filename);
end


%uniquify
E = E(E(:, 1) < E(:, 2), :);

figure
hold on
line([V(E(:,1),1)';V(E(:,2),1)'],[V(E(:,1),2)';V(E(:,2),2)'], 'Color', [0 0 1]);
hold off

levelManager = leveller(V, E);

for i = 1:levelManager.maxLevel
    Einds = levelManager.getEdgesBelowLevel(i);
    levelE = E(Einds, :);
    levelV = levelManager.getVerticesBelowLevel(i);
    
    figure
    hold on
    line([V(levelE(:,1),1)';V(levelE(:,2),1)'],[V(levelE(:,1),2)';V(levelE(:,2),2)'], 'Color', [0 0 1]);
    title(["Input at level " i]);
    hold off
    
    tensions = topologyOptimization(V, levelE, levelV, i);
end

end