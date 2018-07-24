function main(png_filename)
V = [];
E = [];

if (nargin < 1)    
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
    %TO-DO: sometimes this produces a non-stable mesh (one with pointy
    %grounds)
    [Vraw, Fraw] = bwmesh(png_filename);
    [F, V] = reducepatch(Fraw, Vraw, 25);
    
    V = [V(:,1) V(:,2)];
    E = boundary_faces(F);
     
    figure
    tsurf(F, V);
    
    %TO-DO: In order for inpolygon to work, we need to pass
    %in the edges as a polygon, ie. ordered set of edges.
    %So we need to re-order them such that they form something like
    %[a b; b c; c d; ... ]
end

figure
hold on
line([V(E(:,1),1)';V(E(:,2),1)'],[V(E(:,1),2)';V(E(:,2),2)'], 'Color', [0 0 1]);
hold off

levelManager = leveller(V, E);

for i = 1:levelManager.maxLevel
    Einds = levelManager.getEdgesBelowLevel(i);
    levelE = E(Einds, :);
    levelV = levelManager.getVerticesBelowLevel(i);
    
    tensions = topologyOptimization(V, E, levelE, levelV, i);
end

end