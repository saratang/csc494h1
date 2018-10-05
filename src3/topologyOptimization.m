%INPUT
% V - Vertices stored as nx3 matrix
% E - edges stores as mx2 matrix
% levelV - indices of vertices that are in the current level
% levelNum is the current level. I guess I'll get rid of it later since I just use
% it to label figures.
function  dataOpt = topologyOptimization(V,E,levelV, levelNum, F)

    function B = constitutiveModel(vol, ym, eVec)
        %vol is a vector of truss volumes
        %ym young's modulus
        %eVec array of edge vectors
        
        %B matrix (constitutive model) -- this is block diagonal
        iRows = reshape([1:3*numEdges;1:3*numEdges;1:3*numEdges], 9*numEdges,1);
        iCols = reshape([1:3:3*numEdges;
                         2:3:3*numEdges;
                         3:3:3*numEdges;
                         1:3:3*numEdges;
                         2:3:3*numEdges;
                         3:3:3*numEdges;
                         1:3:3*numEdges;
                         2:3:3*numEdges;
                         3:3:3*numEdges;], 9*numEdges,1);
        
        %ordering of constitutive model values [11 12 13 21 22 23 31 32 33]
        lk2 = sum((eVec(1:end,:).^2),2);
        x = (ym./lk2).*(vol.^2);
        vals = reshape([(eVec(:,1).*eVec(:,1))';
                        (eVec(:,1).*eVec(:,2))';
                        (eVec(:,1).*eVec(:,3))';
                        (eVec(:,2).*eVec(:,1))';
                        (eVec(:,2).*eVec(:,2))';
                        (eVec(:,2).*eVec(:,3))';
                        (eVec(:,3).*eVec(:,1))';
                        (eVec(:,3).*eVec(:,2))';
                        (eVec(:,3).*eVec(:,3))'], 9*numEdges, 1);
        vals = reshape(repmat(x', 9, 1), 9*numEdges, 1).*vals;
        
        B = sparse(iRows, iCols, vals);
        
    end

    function WE = reindex_edges(reindex_map, levelV, E)
        % re-index the edges according to the index of each vertex in
        % levelV instead of the value at the index
        
        % basically I don't want to use the find function all the time
        % so I'm going to make a map of Vind->ind in levelV
        
        WE = E;
        
        for vid = 1:numel(E)
            if isKey(reindex_map, E(vid))
                WE(vid) = reindex_map(E(vid));
            else
                reindex_map(E(vid)) = find(levelV == E(vid));
                WE(vid) = reindex_map(E(vid));
            end
        end
    end

adjacencyMatrix = zeros(size(V,1), size(V,1));
adjacencyMatrix(sub2ind(size(adjacencyMatrix),E(:,1), E(:,2))) = 1;
adjacencyMatrix(sub2ind(size(adjacencyMatrix),E(:,2), E(:,1))) = 1;

adjacencyMatrix = adjacencyMatrix - eye(size(V,1));

newTrusses = adjacencyMatrix;
K = [];
At = [];
valid = 0;
Eorig = E;
Enew = E;

while ~valid
    %vertex variables are flattened as [v1x, v1y; v2x, v2y ......]
    Evec = V(E(:,2),:) - V(E(:,1),:);
    
    %Build K matrix
    numEdges = size(E,1);
    numVerts = size(V,1);
    numVars = 3*numVerts; %x,y and z position for each truss
    
    %A matrix (accumulation, distribution)
    % i i+1 i+2 i i+1 i+2
    % j j+1 j+2 j+3 j+4 j+5
    % 1   1   1  -1  -1  -1
    val = repmat([-1 1 -1 1 -1 1], 1, numEdges);
    
    iRow = reshape(repmat(1:(3*numEdges), 2,1),6*numEdges,1);
    E1ind = [3*E(:,1)-2 3*E(:,1)-1 3*E(:,1)]';
    E2ind = [3*E(:,2)-2 3*E(:,2)-1 3*E(:,2)]';
    iCol = reshape([E1ind(:) E2ind(:)]',6*numEdges,1);
    At = sparse(iRow, iCol, val, 3*numEdges, 3*numVerts);
        
    B = constitutiveModel(ones(size(E,1), 1), 1, Evec);
    K = At'*B*At;
    
    floorVerts = find(V(:,3) == min(V(:,3)));
    isolatedVerts = setdiff(1:size(V, 1), levelV).';
    numConstraints = numel(floorVerts);
    P = eye(numel(V), numel(V));
    P([3*[floorVerts; isolatedVerts]-2;
       3*[floorVerts; isolatedVerts]-1;
       3*[floorVerts; isolatedVerts]], :) = [];
    Kcon = P*K*P';
    
    % Get minimal absolute eigenvalue.
    min_eig = min(abs(eig(Kcon)));
    disp(min_eig);
    valid = (min_eig > sqrt(eps));
    
    if ~valid
        E = Eorig;
        newTrusses = newTrusses * adjacencyMatrix;

        %extract new edges
        for ii=1:size(V,1)
            temp = find(newTrusses(ii,:) ~=0);

            for jj = 1:numel(temp)
                if ii ~= temp(jj)   
                    Enew = [Enew; ii temp(jj)];
                end
            end
        end

        Enew = unique([Enew; Enew(:,2) Enew(:,1)], 'rows');
        %add some new edges to make the structure statically over determined

        %strip non unique edges
        [Ecommon1, IA, IB] = intersect(E, Enew, 'rows');
        Enew(IB,:) = [];

        [Ecommon1, IA, IB] = intersect(E, [Enew(:,2) Enew(:,1)], 'rows');
        Enew(IB,:) = [];
        E = [E; Enew];
    end

end

reindex_map = containers.Map('Keytype', 'double', 'Valuetype', 'double');
WV = V(levelV, :, :);
WEorig = reindex_edges(reindex_map, levelV, Eorig);
WEnew = reindex_edges(reindex_map, levelV, Enew);

[Vorig, Forig] = wire_mesh(WV, WEorig, 'Thickness', 0.05);
[Vnew, Fnew] = wire_mesh(WV, WEnew, 'Thickness', 0.05);

figure
ax1 = subplot(1,4,[1 2]);
hold on

trimesh(Fnew, Vnew(:,1,:), Vnew(:,2,:), Vnew(:,3,:), 'edgecolor', 'r');
trimesh(Forig, Vorig(:,1,:), Vorig(:,2,:), Vorig(:,3,:), 'edgecolor', 'b');

title(ax1, ["Overconstrained input at level " levelNum]);
axis equal
view(3);
hold off


fExt = zeros(size(V,1),3);
%fExt(forceVerts, 2) = -0.1;

fExt = reshape(transpose(fExt), numVars, 1);

%The "downwards" direction. Which in this case is the third coordinate.
fExt(3:3:end) = -9.8;

% NOISE
%fExt(1:2:end) = (rand(size(fExt, 1) / 2, 1) - 0.5) * 2;

%Project out fixed constraints from matrices
%detect nodes on "floor"
%for floors
floorVerts = find(V(:,3) == min(V(:,3)));

% isolated vertices ie. vertices that are not attached to edges in this
% level
isolatedVerts = setdiff(1:size(V, 1), levelV).';
numConstraints = numel(floorVerts);
P = eye(numel(V), numel(V));
P([3*[floorVerts; isolatedVerts]-2;3*[floorVerts; isolatedVerts]-1;3*[floorVerts; isolatedVerts]], :) = [];

%define J = A'
J = At;

%remove them from all systems
Kcon = P*K*P';
Jcon = J*P';
fcon = P*fExt;

%Truss direction matrix
bIndices = repmat(1:size(E,1), 3,1);
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

% The cost is the weighted sum of all tensions + compressions
% f is a horizontal vector such that:
%   - roof edges have weight 0... because we always want to include them
%   - overconstrained edges have weight 1.5 or 2 or something... because
%     we want to avoid them if possible
% CAVEAT 1: assign the roof weights first so that we can overwrite them
% with the overconstrained ones if necessary.
f = ones(size(Aeq, 2), 1);

% roofVerts. We want roofVerts because we want to assign roof edges a
% super-low cost. For now let's do it for the "absolute roof", but we can
% change this later to the incremental roof.
roofVerts = find(V(:,3) == max(V(:,3)));
[X, Y] = meshgrid(roofVerts, roofVerts);
roofEdges = intersect([X(:) Y(:)], E, 'rows');

% TODO: ask Dave what the better way to do this is lol
% TODO: how do I do runtime analysis in Matlab
roofInds = find(ismember(E, roofEdges, 'rows'));
roofInds = [roofInds roofInds + size(E,1)];
f(roofInds) = 0;

% TODO: We want to distinguish between overconstrainedInds that are inside
% our structure and those outside.
overconstrainedInds = [size(E, 1) - size(Enew, 1) + 1:size(E, 1)];
overconstrainedInds = [overconstrainedInds overconstrainedInds + size(E,1)];
f(overconstrainedInds) = 1e3;

% TODO: Find edges that are not in the polygon.
% For now, assume we have an outline and not a triangle mesh.
% Do this by querying the midpoint of each edge and seeing if it's in the
% polygon.

% Use winding_number from gp-toolbox. To do this, we need a triangle mesh.


% FORMAT: winding_number(V, F, O) where
%  - V is a 3-list of vertex poistions
%  - F is a 3-list of triangle indices
%  - O is a 3-list of queries
%  polygon

midpoints = ((V(E(:,1),:) + V(E(:,2),:)) / 2);
onequarter = ((3 * V(E(:,1),:) + V(E(:,2),:)) / 4);
threequarter = ((V(E(:,1),:) + 3 * V(E(:,2),:)) / 4);

in1 = winding_number(V, F, midpoints);
in2 = winding_number(V, F, onequarter);
in3 = winding_number(V, F, threequarter);

% Now we select the indices of those that are outside and set them to an
% insanely high weight
outsideEdgeInds = find(in1 < 1e-04 | in2 < 1e-04 | in3 < 1e-04);
outsideEdgeInds = [outsideEdgeInds outsideEdgeInds + size(E,1)];
f(outsideEdgeInds) = 1e6;

% Also add the constraint that all tensions, compressions >= 0
A = eye(size(f, 1));
b = zeros(size(f, 1), 1);

%X should be a col. vector with top half tension, lower half compressions
X = linprog(f, -A, b, Aeq, Beq);
%disp(X);

%We want to draw tension - compression
%waitwaitwait I'm seeing compressions. Also are my dimensions off?
tensions = X(1:size(X, 1) / 2);
compressions = X(size(X, 1) / 2 + 1 : end);
%disp(tensions);

stress = abs(tensions - compressions);

%HACK
threshold = 1e-04;

dataOpt = stress;

%plot result
%I imagine you would map the tensions (top half of X)
%back to the original "edges"
%don't plot the zero ones
%plot the others with different colours corresponding to different tensions

disp("stress");
disp(stress);

stressedEdges = E(stress > threshold,:);

WEfinal = reindex_edges(reindex_map, levelV, [stressedEdges; Eorig]);
[Vfinal, Ffinal] = wire_mesh(WV, WEfinal, 'Thickness', 0.05);

ax2 = subplot(1,4,[3 4]);
hold on
trimesh(Ffinal, Vfinal(:,1,:), Vfinal(:,2,:), Vfinal(:,3,:), 'edgecolor', 'b');
title(ax2,["Optimized at level " levelNum]);
axis equal;
view(3)
hold off

end

