

%INPUT
% V - Vertices stored as nx3 matrix
% E - edges stores as mx2 matrix
% levelV - indices of vertices that are in the current level
% levelNum is the current level. I guess I'll get rid of it later since I just use
% it to label figures.
function  dataOpt = topologyOptimization(V,allE,levelE,levelV, levelNum)

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

% TODO: This is not generating the trusses I want. How can I generate the
% trusses I want?
adjacencyMatrix = zeros(size(V,1), size(V,1));
adjacencyMatrix(sub2ind(size(adjacencyMatrix),levelE(:,1), levelE(:,2))) = 1;
adjacencyMatrix(sub2ind(size(adjacencyMatrix),levelE(:,2), levelE(:,1))) = 1;


adjacencyMatrix = adjacencyMatrix - eye(size(V,1));

newTrusses = adjacencyMatrix;
K = [];
At = [];
valid = 0;

while ~valid
    newTrusses = newTrusses * adjacencyMatrix;
    
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
    [Ecommon1, IA, IB] = intersect(levelE, Enew, 'rows');
    Enew(IB,:) = [];
    
    [Ecommon1, IA, IB] = intersect(levelE, [Enew(:,2) Enew(:,1)], 'rows');
    Enew(IB,:) = [];
    Eorig = levelE;
    levelE = [levelE; Enew];
    
    %vertex variables are flattened as [v1x, v1y; v2x, v2y ......]
    Evec = V(levelE(:,2),:) - V(levelE(:,1),:);
    
    %Build K matrix
    numEdges = size(levelE,1);
    numVerts = size(V,1);
    numVars = 2*numVerts; %x and y position for each truss
    numVarsUncompressed = 2*numEdges;
    
    %A matrix (accumulation, distribution)
    % i i+1 i+2 i i+1 i+2
    % j j+1 j+2 j+3 j+4 j+5
    % 1   1   1  -1  -1  -1
    val = repmat([-1 1 -1 1], 1, numEdges);
    iRow = reshape(repmat(1:numVarsUncompressed, 2,1),2*numVarsUncompressed,1);
    E1ind = [2*levelE(:,1)-1 2*levelE(:,1)]';
    E2ind = [2*levelE(:,2)-1 2*levelE(:,2)]';
    iCol = reshape([E1ind(:) E2ind(:)]',2*numVarsUncompressed,1);
    At = sparse(iRow, iCol, val, 2*numEdges, numVars);
    
    forceVerts = (find(and(abs(V(:,1) - max(V(:,1)))<1e-8,V(:,2) == 0)));
    
    B = constitutiveModel(ones(size(levelE,1), 1), 1, Evec);
    K = At'*B*At;
    
    floorVerts = find(V(:,2) == min(V(:,2)));
    isolatedVerts = setdiff(1:size(V, 1), levelV).';
    numConstraints = numel(floorVerts);
    P = eye(numel(V), numel(V));
    P([2*[floorVerts; isolatedVerts]-1;2*[floorVerts; isolatedVerts]], :) = [];
    Kcon = P*K*P';
    
    % CONSULT DAVE
    % Get minimal absolute eigenvalue.
    min_eig = min(abs(eig(Kcon)));
    disp(min_eig);
    valid = (min_eig > sqrt(eps));
    
    if ~valid
        levelE = Eorig;
    end
end
    
figure
ax1 = subplot(1,4,[1 2]);
hold on
line(ax1, [V(Enew(:,1),1)';V(Enew(:,2),1)'],[V(Enew(:,1),2)';V(Enew(:,2),2)'], 'Color', [1 0 0]);
line(ax1, [V(Eorig(:,1),1)';V(Eorig(:,2),1)'],[V(Eorig(:,1),2)';V(Eorig(:,2),2)'], 'Color', [0 0 1]);
title(ax1, ["Overconstrained input at level " levelNum]);
axis square;
hold off

fExt = zeros(size(V,1),2);
%fExt(forceVerts, 2) = -0.1;

fExt = reshape(transpose(fExt), numVars, 1);
fExt(2:2:end) = -9.8;

% NOISE
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

%define J = A'
J = At;

%remove them from all systems
Kcon = P*K*P';
Jcon = J*P';
fcon = P*fExt;

%Truss direction matrix
bIndices = repmat(1:size(levelE,1), 2,1);
bI = 1:numel(Evec);
bJ = bIndices(:);

norm_E = sqrt(sum(Evec.^ 2, 2));     %// Calculate Euclidean length
norm_E(norm_E < eps) = 1;            %// Avoid division by zero
dE = bsxfun(@rdivide, Evec, norm_E); %// Normalize
T = sparse(bI, bJ, reshape(transpose(dE), numel(dE),1));

%Build Dual problem matrices
Jcon = T'*Jcon;

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
roofVerts = find(V(:,2) == max(V(:,2)));
[X, Y] = meshgrid(roofVerts, roofVerts);
roofEdges = intersect([X(:) Y(:)], levelE, 'rows');

% TODO: ask Dave what the better way to do this is lol
% TODO: how do I do runtime analysis in Matlab
roofInds = find(ismember(levelE, roofEdges, 'rows'));
roofInds = [roofInds roofInds + size(levelE,1)];
f(roofInds) = 0;

overconstrainedInds = [size(levelE, 1) - size(Enew, 1) + 1:size(levelE, 1)];
overconstrainedInds = [overconstrainedInds overconstrainedInds + size(levelE,1)];
f(overconstrainedInds) = 1e3;

% TODO: Find edges that are not in the polygon.
% For now, assume we have an outline and not a triangle mesh.
% Do this by querying the midpoint of each edge and seeing if it's in the
% polygon.

% FORMAT: inpolygon(xq, yq, xv, yv) where
%  - xq, yq are the coordinates we want to query
%  - xv, yv are row vectors such that (xvi, yvi) defines a vertex of the
%  polygon

polygonV = polygonize(V, allE);
xv = polygonV(:,1)';
yv = polygonV(:,2)';
midpoints = ((V(levelE(:,1),:) + V(levelE(:,2),:)) / 2)';
onequarter = ((3 * V(levelE(:,1),:) + V(levelE(:,2),:)) / 4)';
threequarter = ((V(levelE(:,1),:) + 3 * V(levelE(:,2),:)) / 4)';
in1 = inpolygon(midpoints(1,:), midpoints(2,:), xv, yv);
in2 = inpolygon(onequarter(1,:), onequarter(2,:), xv, yv);
in3 = inpolygon(threequarter(1,:), threequarter(2,:), xv, yv);

% Now we select the indices of those that are outside and set them to an
% insanely high weight
outsideEdgeInds = find(~in1 | ~in2 | ~in3);
outsideEdgeInds = [outsideEdgeInds outsideEdgeInds + size(levelE,1)];
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

stressedEdges = levelE(stress > threshold,:);
edgeSrcs = [stressedEdges(:,1); Eorig(:,1)];
edgeDsts = [stressedEdges(:,2); Eorig(:,2)];

ax2 = subplot(1,4,[3 4]);
hold on
line(ax2, ...
     [V(edgeSrcs,1)'; ...
      V(edgeDsts,1)'], ...
     [V(edgeSrcs,2)'; ...
      V(edgeDsts,2)'], ...
     'Color', [0 0 1]);
title(ax2,["Optimized at level " levelNum]);
% WARNING: specific to the current design. But I want the axes to match
% so...
axis square;
hold off

end

