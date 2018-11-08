function main(obj_filename)
V = [];
E = [];
F = [];

if (nargin < 1)
    %IMMA DRAW THE T. In 3d!
    V = [0 0 2;
         0 0 3;
         1 0 3;
         2 0 3;
         3 0 3;
         3 0 2;
         2 0 2;
         2 0 1;
         2 0 0;
         1 0 0;
         1 0 1;
         1 0 2;
         0 1 2;
         0 1 3;
         1 1 3;
         2 1 3;
         3 1 3;
         3 1 2;
         2 1 2;
         2 1 1;
         2 1 0;
         1 1 0;
         1 1 1;
         1 1 2];
          
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
         12 1;
         13 14;
         14 15;
         15 16;
         16 17;
         17 18;
         18 19;
         19 20;
         20 21;
         21 22;
         22 23;
         23 24;
         24 13;
         1 13;
         2 14;
         5 17;
         6 18;
         7 19;
         9 21;
         10 22;
         12 24;];
     
%      disp('Contents of workspace before loading file:');
%      whos
%      disp('Contents of results_sigAsia_1_4.mat:');
%      whos('-file', '../meshes/results_sigAsia_1_4.mat');
    
else % expect and OBJ from Blender
    [V, F] = readOBJ(obj_filename);
    % for some reason the display is rotated so the first thing
    % I'm going to do is un-rotate it
    rot = [1 0 0; 0 0 1; 0 -1 0];
    V = V * rot;
    
    figure
    hold on
    tsurf(F, V);
    axis equal;
    view(3);
    hold off
    
    E = edges(F);
end

figure
hold on
line([V(E(:,1),1)'; ...
      V(E(:,2),1)'],...
     [V(E(:,1),2)'; ...
      V(E(:,2),2)'],...
     [V(E(:,1),3)'; ...
      V(E(:,2),3)']);
axis equal;
view(3);
hold off

% figure
% hold on
% line([V(E(:,1),1)';V(E(:,2),1)'],[V(E(:,1),2)';V(E(:,2),2)'], 'Color', [0 0 1]);
% hold off

levelManager = leveller(V, E);

colors = [ 0 0 1; 0 1 0; 1 0 0];

figure
for i = 0:levelManager.maxLevel
    Einds = levelManager.getEdgesInLevel(i);
    levelE = E(Einds, :);
    
    hold on
    line([V(levelE(:,1),1)';V(levelE(:,2),1)'],[V(levelE(:,1),2)';V(levelE(:,2),2)'], [V(levelE(:,1),3)'; V(levelE(:,2),3)'], 'Color', colors(mod(i, 3) + 1, :));
    axis equal;
    view(3);
    hold off
    
end

for i = 1:levelManager.maxLevel
    Einds = levelManager.getEdgesBelowLevel(i);
    levelE = E(Einds, :);
    levelV = levelManager.getVerticesBelowLevel(i);
    
     figure
     hold on
     line([V(levelE(:,1),1)';V(levelE(:,2),1)'],[V(levelE(:,1),2)';V(levelE(:,2),2)'], [V(levelE(:,1),3)'; V(levelE(:,2),3)'], 'Color', [0 0 1]);
     axis equal;
     view(3);
     hold off
    
     %tensions = topologyOptimization(V, levelE, levelV, i, F);
end

end