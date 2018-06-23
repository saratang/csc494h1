classdef leveller
    properties        
        % This is a container of level-to-array-of-edges
        levelsToEdges
        
        % This is a container of level-to-array-of-vertices
        levelsToVertices
        
        % This is a map of edge-to-level
        edgesToLevels
        
        % This is a map of vertex-to-level
        verticesToLevels
        
        % Number representing top-level
        maxLevel
    end
    methods
        % Given the mesh and the ground nodes, populate the levels using BFS
        function obj = leveller(vertices, edges)
            groundVertices = find(vertices(:,2) == min(vertices(:,2)));
            disp(groundVertices);
            
            obj.levelsToEdges = containers.Map('KeyType', 'double', 'ValueType', 'any');
            obj.levelsToVertices = containers.Map('KeyType', 'double', 'ValueType', 'any');
            obj.maxLevel = 0;
            
            function ind = getEdgeIndex(edges, edge)
                ind = find((edges(:, 1) == edge(1) & edges(:, 2) == edge(2)) | (edges(:, 1) == edge(2) & edges(:, 2) == edge(1)));
            end
            
            [X, Y] = meshgrid(groundVertices, groundVertices);
            %cartesian product of ground vertices intersect edges
            groundEdges = intersect([X(:) Y(:)], edges, 'rows');
            obj.edgesToLevels = Inf(1, size(edges, 1));
            for row = 1:size(groundEdges, 1)
                edge = groundEdges(row, :);
                edgeID = getEdgeIndex(edges, edge);
                obj.edgesToLevels(edgeID) = 0;
            end
            
            obj.verticesToLevels = Inf(1, size(vertices, 1));
            for row = 1:size(groundVertices, 1)
                vID = groundVertices(row);
                obj.verticesToLevels(vID) = 0;
            end
            
            G = graph(table(edges, 'variableNames', {'EndNodes'}));
            for i = 1:size(groundVertices, 1)
                groundV = groundVertices(i);
                T = bfsearch(G, groundV, {'edgetonew', 'edgetodiscovered'});
                curLevel = 1;
                updateLevel = NaN;
                for eid = 1:size(T, 1)
                    event = T.Event(eid);
                    edge = T.Edge(eid, :);
                    node = T.Node(eid);
                    if (edge(1) == updateLevel)
                        updateLevel = NaN;
                        curLevel = curLevel + 1;
                    end
                    
                    if event == 'edgetonew'
                        if (isnan(updateLevel))
                            updateLevel = edge(2);
                        end
                        edgeID = getEdgeIndex(edges, edge);
                        obj.edgesToLevels(edgeID) = min(obj.edgesToLevels(edgeID), curLevel);
                        obj.verticesToLevels(edge(2)) = min(obj.verticesToLevels(edge(2)), curLevel);
                    end
                    if event == 'edgetodiscovered'
                        edgeID = getEdgeIndex(edges, edge);
                        obj.edgesToLevels(edgeID) = min(obj.edgesToLevels(edgeID), curLevel);
                    end
                end
                disp(T);
            end
            
            % invert the map to get {level : array of edges}
            % MAY NOT NEED THIS MAP
            for eid = 1:size(obj.edgesToLevels, 2)
                level = obj.edgesToLevels(eid);
                
                if (level > obj.maxLevel)
                    obj.maxLevel = level;
                end
                
                if (isKey(obj.levelsToEdges, level))
                    obj.levelsToEdges(level) = [obj.levelsToEdges(level); eid];
                else
                    obj.levelsToEdges(level) = eid;
                end
            end
            
            for vid = 1:size(obj.verticesToLevels, 2)
                level = obj.verticesToLevels(vid);
                
                fprintf("%d: lvl %d\n", vid, level);
                
                if (isKey(obj.levelsToVertices, level))
                    obj.levelsToVertices(level) = [obj.levelsToVertices(level); vid];
                else
                    obj.levelsToVertices(level) = vid;
                end
            end
        end
        
        % Return the set of edge indices at or below levelNum
        function E = getEdgesBelowLevel(obj, levelNum)
            assert(levelNum <= obj.maxLevel, "level number %d is larger than the max level, %d", levelNum, obj.maxLevel);
            
            E = [];
            for level = 0:levelNum
                E = [E; obj.levelsToEdges(level)];
            end
        end
        
        % Return a list of vertex indices at or below levelNum
        function V = getVerticesBelowLevel(obj, levelNum)
            assert(levelNum <= obj.maxLevel, "level number %d is larger than the max level, %d", levelNum, obj.maxLevel);
            
            V = [];
            for level = 0:levelNum
                if (isKey(obj.levelsToVertices, level)) % It is possible for a level to exist with no "new" nodes
                    V = [V; obj.levelsToVertices(level)];
                end
            end
        end
        
        % Return the set of bars in levelNum
        function E = getEdgesInLevel(obj, levelNum)
            assert(levelNum <= obj.maxLevel, "level number %d is larger than the max level, %d", levelNum, obj.maxLevel);
            
            E = obj.levelsToEdges(levelNum);
        end
        
        function V = getVerticesInLevel(obj, levelNum)
            assert(levelNum <= obj.maxLevel, "level numer %d is larger than max level, %d", levelNum, obj.maxLevel);
            
            V = obj.levelsToVertices(levelNum);
        end
    end
end