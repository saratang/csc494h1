function  vOut = trussOpt(data)

%optimize a sparse struss structure
if nargin < 1
    data = [];
    data.V = [];
    data.fExt = []; 
    data.Curves = [];
end

%step one initial mesh
if isempty(data.V)
    %sample data
    nX = 10;
    
    % make a coordinate grid of size 10x10
    [X,Y] = meshgrid(linspace(0,1,nX), linspace(0,1,nX));
    V = [X(:), Y(:)];
    
    % triangulate the grid. It's 1-based indexing btw.
    F = delaunayTriangulation(V(:,1), V(:,2));
    
    % get the list of edges (edge = vertex-pair)
    E = edges(F);
    
    % probably gets it in a nicer form
    Evec = V(E(:,2),:) - V(E(:,1),:); 
else
    %use mesh from input data
    V = data.Vopt;
    F = data.Fopt;
    E = edges(F);
    Evec = V(E(:,2),:) - V(E(:,1),:); 
    
    %nX = sum(V(:,1) == max(V(:,1)));
    fExt = 0*V;
    nX = sum(V(:,1) == max(V(:,1)));
    cornerIndex = and(V(:,1)==max(V(:,1)),V(:,2)==min(V(:,2)));
    fExt(find(cornerIndex==1)+int32(nX/2),:) = [1, -10];
    %fExt((nX-1)*nX+1 + int32(nX/2), 2) = -10;
    %fExt((nX-1)*nX+1 + int32(nX/2), 1) = 1;

    %fix nodes on left hand side, [1/2,2/3] of the way up
    iFloor = (V(:,1) == min(V(:,1)));
    
    %use simple boundary conditions or read from data in this case
    %if isempty(data.fExt)
        %apply forces on top
        %fExt = 0*V;
        %fExt(V(:,2) == max(V(:,2)), 2) = -20;
    %else
        %fExt = data.fExt;
    %end
    
    %if isempty(data.fixed)
        %fix points on the bottom
        %iFloor = (V(:,2) == min(V(:,2)));
    %else
     %   iFloor = data.fixed;
    %end
    
    if(isempty(data.Curves))
        
    else
        curves = data.Curves; 
    end
end

numEdges = size(E,1);
numVerts = size(V,1); 
numVarsUncompressed = 2*numEdges;
numVars = 2*numVerts;

%A matrix (accumulation, distribution)
% i i+1 i+2 i i+1 i+2 
% j j+1 j+2 j+3 j+4 j+5 
% 1   1   1  -1  -1  -1 
val = repmat([-1 1 -1 1], 1, numEdges); 
iRow = reshape(repmat(1:numVarsUncompressed, 2,1),2*numVarsUncompressed,1);

%The apostrophe is complex conjugate transpose.
E1ind = [2*E(:,1)-1 2*E(:,1)]';
E2ind = [2*E(:,2)-1 2*E(:,2)]';
iCol = reshape([E1ind(:) E2ind(:)]',2*numVarsUncompressed,1);
At = sparse(iRow, iCol, val, 2*numEdges, numVars); 

    function B = constitutiveModel(vol, ym, eVec)
        %vol is a vector of truss volumes
        %ym young's modulus
        %eVecarray of edge vectors
        
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

B = constitutiveModel(ones(size(E,1), 1), 1000, Evec);
K = At'*B*At;

if isempty(data.V) 
    %use parameteric test cases for this 
    %simple test case
    %fExt (external forces)
    fExt = 0*V;
    %fExt(90, 1) = 9.8; 
    %fExt(90, 2) = -9.8; 
    fExt(V(:,2) == max(V(:,2)),2) = -9.8;

    %fix floor nodes
    iFloor = (V(:,2) == min(V(:,2))); 
  

    %hanging design
    %add downward load at center of right most edge
    fExt = 0*V;
    fExt((nX-1)*nX+1 + int32(nX/2), 2) = -10;
    fExt((nX-1)*nX+1 + int32(nX/2), 1) = 1;

    %fix nodes on left hand side, [1/2,2/3] of the way up
    minY = min(V(:,2));
    maxY = max(V(:,2));
    lower = 0.0*(maxY - minY) + minY;
    upper = 1.0*(maxY - minY) + minY;
    iFloor = and(and((V(:,2) > lower), (V(:,2)<upper)), V(:,1) == min(V(:,1)));

    %bridge design
else
    %get fixed nodes, forces from data
end
fExt = reshape(transpose(fExt), numel(fExt), 1);
iFloor = [iFloor'; iFloor'];
iFloor = iFloor(:);
fExt(iFloor) = [];

%plot results 
%h = figure;

%updated vertices
%Vnew = V + uComplete; 
%line([Vnew(E(:,1),1)'; Vnew(E(:,2),1)'], [Vnew(E(:,1),2)'; Vnew(E(:,2),2)'], 'Color', 'blue');

%optimization stuff

    %cost function which is the negative compliance of the structure
    function [c, g] = compliance(vol, fExt, ym, eVec, At, fixed)
         numVerts = numel(fExt);
         B = constitutiveModel(vol, ym, eVec);
         K = At'*B*At;
         K(fixed,:) = [];
         K(:, fixed) = [];
         u = K\fExt;
         uFilled = zeros(numVerts,1);
         uFilled(~fixed) = u;
         c = 0.5.*fExt'*u;
         u = uFilled;
         
         %implement gradient
         if nargin > 1
            du = At*u; 
            volNew = (2.*vol).^(1/2);
            dK = constitutiveModel(volNew, ym, eVec);
            parfor i=1:size(eVec, 1)
                g(i) = -0.5.*du((2*(i-1)+1):(2*(i-1)+2))'*dK((2*(i-1)+1):(2*(i-1)+2),(2*(i-1)+1):(2*(i-1)+2))*du((2*(i-1)+1):(2*(i-1)+2));
            end
            g = g';
         end
    end

    vTotal = 1;
    
    %constraints, volume >= 0, volume < max volume allowed for a truss
    lb = 0.001.*ones(size(E,1), 1);
    %ub = (size(E,1)/2)*ones(size(E,1),1);
    ub = [];
    %total amount of volume allowed
    Aeq = ones(1, size(E,1));
    beq = size(E,1);
    
    %initial Volume 
    v0 = (vTotal./size(E,1)).*ones(size(E,1),1);

    %trusses that are part of the user design must have some amount of volume
    
    % This is the (cheat) algorithm that computes the min-cost truss
    options = optimoptions(@fmincon, 'Algorithm', 'interior-point', 'Hessian', 'lbfgs', 'MaxFunEvals', 1000000, 'MaxIter', 2000, 'Display', 'iter', 'GradObj','on', 'DerivativeCheck', 'off', 'InitBarrierParam', 100);
    vDone = fmincon(@(x) compliance(x, fExt,1000, Evec, At, iFloor), v0, [], [], Aeq, beq,lb, ub, [], options); 
    vOut = vDone;
    vDone(vDone < 0.01*max(vOut)) = 0;
    vDone(vDone > 0) = 1;
    %draw trusses but color with vDone
    Vnew = V; 
    figure
    hold on
    
  
    for ii=1:size(E,1)
        alpha = (vDone(ii)./max(vDone));
        line([Vnew(E(ii,1),1)'; Vnew(E(ii,2),1)'], [Vnew(E(ii,1),2)'; Vnew(E(ii,2),2)'], 'Color', [0.0,0.0,1.0].*alpha+[1.0,1.0,1.0].*(1-alpha), 'LineWidth', 4.0.*vOut(ii)./max(vOut));
    end
    
    plot(V(iFloor(1:2:end),1), V(iFloor(1:2:end),2), 'r*');
    hold off

end

