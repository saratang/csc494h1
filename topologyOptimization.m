%INPUT
% V - Vertices stored as nx3 matrix
% E - edges stores as mx2 matrix
% levelV - indices of vertices that are in the current level
% i is the current level. I guess I'll get rid of it later since I just use
% it to label figures.
function  dataOpt = topologyOptimization(V,E,levelV, levelNum)

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

figure
hold on
line([V(E(:,1),1)';V(E(:,2),1)'],[V(E(:,1),2)';V(E(:,2),2)'], 'Color', [0 0 1]);
title(["Overconstrained input at level " levelNum])
hold off

%vertex variables are flattened as [v1x, v1y; v2x, v2y ......]
Evec = V(E(:,2),:) - V(E(:,1),:);

%Build K matrix
numEdges = size(E,1);
numVerts = size(V,1);
numVars = 2*numVerts; %x and y position for each truss
numVarsUncompressed = 2*numEdges;

%A matrix (accumulation, distribution)
% i i+1 i+2 i i+1 i+2
% j j+1 j+2 j+3 j+4 j+5
% 1   1   1  -1  -1  -1
val = repmat([-1 1 -1 1], 1, numEdges);
iRow = reshape(repmat(1:numVarsUncompressed, 2,1),2*numVarsUncompressed,1);
E1ind = [2*E(:,1)-1 2*E(:,1)]';
E2ind = [2*E(:,2)-1 2*E(:,2)]';
iCol = reshape([E1ind(:) E2ind(:)]',2*numVarsUncompressed,1);
At = sparse(iRow, iCol, val, 2*numEdges, numVars);

    function B = constitutiveModel(vol, ym, eVec)
        %vol is a vector of truss volumes
        %ym young's modulus
        %eVec array of edge vectors
        
        %B matrix (constitutive model) -- this is block diagonal
        iRows = reshape([1:2*numEdges;1:2*numEdges], 4*numEdges,1);
        iCols = reshape([1:2:2*numEdges;2:2:2*numEdges;1:2:2*numEdges;2:2:2*numEdges], 4*numEdges,1);
        
        %ordering of constitutive model values [11 12 21 22]
        lk2 = sum((eVec(1:end,:).^2),2);
        x = (ym./lk2).*(vol.^2);
        vals = reshape([(eVec(:,1).*eVec(:,1))'; (eVec(:,1).*eVec(:,2))'; (eVec(:,2).*eVec(:,1))'; (eVec(:,2).*eVec(:,2))'], 4*numEdges, 1);
        vals = reshape(repmat(x', 4, 1), 4*numEdges, 1).*vals;
        
        B = sparse(iRows, iCols, vals);
        
    end

forceVerts = (find(and(abs(V(:,1) - max(V(:,1)))<1e-8,V(:,2) == 0)));

B = constitutiveModel(ones(size(E,1), 1), 1, Evec);
K = At'*B*At;
fExt = zeros(size(V,1),2);
%fExt(forceVerts, 2) = -0.1;

fExt = reshape(transpose(fExt), numVars, 1);
fExt(2:2:end) = -9.8;
%fExt(1:2:end) = (rand(size(fExt, 1) / 2, 1) - 0.5) * 2;

%Project out fixed constraints from matrices
%detect nodes on "floor"
%for floors
floorVerts = find(V(:,2) == min(V(:,2)));

% isolated vertices ie. vertices that are not attached to edges in this
% level
isolatedVerts = setdiff(1:size(V, 1), levelV).';
numConstraints = numel(floorVerts);
P = eye(numel(V), numel(V));
P([2*[floorVerts; isolatedVerts]-1;2*[floorVerts; isolatedVerts]], :) = [];

% roofVerts. We want roofVerts because we want to assign roof edges a
% super-low cost. For now let's do it for the "absolute roof", but we can
% change this later to the incremental roof.
roofVerts = find(V(:,2) == max(V(:,2)));
[X, Y] = meshgrid(roofVerts, roofVerts);
roofEdges = intersect([X(:) Y(:)], E, 'rows');

%define J = A'
J = At;

%remove them from all systems
Kcon = P*K*P';
Jcon = J*P';
fcon = P*fExt;

%Truss direction matrix
bIndices = repmat(1:size(E,1), 2,1);
bI = 1:numel(Evec);
bJ = bIndices(:);

norm_E = sqrt(sum(Evec.^ 2, 2));     %// Calculate Euclidean length
norm_E(norm_E < eps) = 1;            %// Avoid division by zero
dE = bsxfun(@rdivide, Evec, norm_E); %// Normalize
T = sparse(bI, bJ, reshape(transpose(dE), numel(dE),1));

%Build Dual problem matrices
Jcon = T'*Jcon;

% %METHOD 1: Solving the linear system.
% %set-up the appropriate coefficient matrix
% coefficients = [Kcon Jcon'; Jcon zeros(size(Jcon,1), size(Jcon,1))];
%
% %set-up the appropriate constant matrix
% constants = [fcon; zeros(size(Jcon,1), size(fcon,2))];
%
% %do the solve
% %UT = coefficients \ constants;
%
% % MAGIC
% %disp(size(coefficients))
% %disp(size(constants))
% UT = quadprog(eye(size(coefficients, 1)), zeros(size(coefficients, 1), 1), [], [], coefficients, constants);
% %U = quadprog(Kcon, -fcon, [], [], Jcon, zeros(size(Jcon,1), size(fcon,2)));
%
% %U = linsolve(coefficients, constants);
%
% %disp(size(UT));
% %disp(UT);
% %disp(U);


% METHOD 2: L1-optimization
% D-matrix looks like this: I | -I

% Given.
Aeq = - (Jcon * inv(Kcon) * Jcon');

% Augment Aeq
D = [eye(size(Aeq, 2)) -eye(size(Aeq, 2))];
Aeq = Aeq * D;

% Given.
Beq = Jcon * inv(Kcon) * fcon;

% The cost is the sum of all tensions + compressions
f = ones(size(Aeq, 2), 1);

% Also add the constraint that all tensions, compressions >= 0
A = eye(size(f, 1));
b = zeros(size(f, 1), 1);

%X should be a col. vector with top half tension, lower half compressions
X = linprog(f, -A, b, Aeq, Beq);
%disp(X);

%waitwaitwait I'm seeing compressions. Also are my dimensions off?
tensions = X(1:size(X, 1) / 2);
%disp(tensions);

dataOpt = tensions;

%plot result
%I imagine you would map the tensions (top half of X)
%back to the original "edges"
%don't plot the zero ones
%plot the others with different colours corresponding to different tensions

disp("V");
disp(V);
disp("E");
disp(E);
disp("tensions");
disp(tensions);

figure
hold on
line([V(E(tensions > 0,1),1)'; ...
    V(E(tensions > 0,2),1)'], ...
    [V(E(tensions > 0,1),2)'; ...
    V(E(tensions > 0,2),2)'], ...
    'Color', [0 0 1]);
title(["Optimized at level " levelNum]);
hold off

end
