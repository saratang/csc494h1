function levelTester
nX = 4;
nY= 4;

[X,Y] = meshgrid(linspace(-5,5,nX), linspace(-2.5,2.5,nY));
V = [X(:), Y(:)];

%build trusses by triangulating nodes
F = delaunayTriangulation(V(:,1), V(:,2));
E = edges(F);

adjacencyMatrix = zeros(size(V,1), size(V,1));
adjacencyMatrix(sub2ind(size(adjacencyMatrix),E(:,1), E(:,2))) = 1;
adjacencyMatrix(sub2ind(size(adjacencyMatrix),E(:,2), E(:,1))) = 1;


adjacencyMatrix = adjacencyMatrix - eye(size(V,1));

newTrusses = adjacencyMatrix*adjacencyMatrix;

Enew = [];

%extract new edges
for ii=1:size(V,1)
    temp = find(newTrusses(ii,:) ~=0);
    
    for jj = 1:numel(temp)
        Enew = [Enew; ii temp(jj)];
    end
end

Enew = unique([Enew; Enew(:,2) Enew(:,1)], 'rows');
%add some new edges to make the structure statically over determined

%strip non unique edges
[Ecommon1, IA, IB] = intersect(E, Enew, 'rows');
Enew(IB,:) = [];

[Ecommon1, IA, IB] = intersect(E, [Enew(:,2) Enew(:,1)], 'rows');
Enew(IB,:) = [];
Eorig = E;
E = [E; Enew];

%uniquify
E = E(E(:, 1) < E(:, 2), :);

myLeveller = leveller(V, E);

for i = 0:myLeveller.maxLevel
    fprintf("Edges in level %d\n", i);
    disp(myLeveller.getEdgesInLevel(i));
end

for i = 0:myLeveller.maxLevel
    fprintf("Edges at or below level %d\n", i);
    disp(myLeveller.getEdgesBelowLevel(i));
end

figure
hold on
for i = 0:myLeveller.maxLevel
    edgeIndices = myLeveller.edgesToLevels == i;
    color = i / myLeveller.maxLevel;
    disp(color);
    line([V(E(edgeIndices,1),1)'; ...
        V(E(edgeIndices,2),1)'], ...
        [V(E(edgeIndices,1),2)'; ...
        V(E(edgeIndices,2),2)'], ...
        'Color', [color color color], ...
        'LineWidth', myLeveller.maxLevel - i + 0.5);
end
hold off
end