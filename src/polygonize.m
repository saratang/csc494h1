% INPUT
% - a mx2 matrix of edge indices
% OUTPUT
% - a re-ordered set of edges of the form
% [a b; b c; c d; ...]
% - a re-ordered set of vertex coordinates
% [xa ya; xb yb; ...]
% Assume each src has a single dest and vice versa
function newVertices = polygonize(V, E)
    edgeDict = containers.Map('KeyType','int32', 'ValueType', 'int32');
    for i = 1:size(E, 1)
        edgeSrc = E(i, 1);
        edgeDst = E(i, 2);
        
        edgeDict(edgeSrc) = edgeDst;
    end
    
    % Start populating newEdges
    newEdges = zeros(size(E));
    newEdges(1,:) = E(1,:);
        
    curSrc = E(1, 1);
    curDst = E(1, 2);
    curInd = 2;
    loopSrc = E(1, 1); % not needed... yet
    remove(edgeDict, curSrc);
    
    % TODO: this won't work for shapes with holes in it
    while ~isempty(edgeDict)
        if ~isKey(edgeDict, curDst)
            curSrc = curDst;
            curDst = loopSrc;
        else
            curSrc = curDst;
            curDst = edgeDict(curDst);
            
            remove(edgeDict, curSrc);
        end
        
        newEdges(curInd,:) = [curSrc curDst];

        curInd = curInd + 1;
    end
    
    % Alternatively, you could reorder everything.
    newVertices = V(newEdges(:, 1),:);
        
end