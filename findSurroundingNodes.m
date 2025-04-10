function surroundingNodes = findSurroundingNodes(boundaryNodesSet, v)
    % Initialise storage structure for all letters
    letters = fieldnames(boundaryNodesSet);
    surroundingNodes = struct();

    % Loop over each letter (C, D, Q)
    for i = 1:length(letters)
        letter = letters{i};
        boundaryNodes = boundaryNodesSet.(letter);
        
        % Initialise an empty array for surrounding nodes
        surroundingNodesForLetter = [];
        
        for k = 1:length(boundaryNodes)
            % Find triangles containing the boundary node
            [row_idx, ~] = find(v == boundaryNodes(k));

            % Get all nodes in those triangles
            nodes_in_triangles = unique(v(row_idx, :));

            % Add new nodes that are not already in boundaryNodes
            surroundingNodesForLetter = unique([surroundingNodesForLetter; setdiff(nodes_in_triangles, boundaryNodes)]);
        end
        
        % Ensure surrounding nodes are within valid range (1 to 256)
        surroundingNodesForLetter = surroundingNodesForLetter(surroundingNodesForLetter >= 1 & surroundingNodesForLetter <= 256);

        % Store surrounding nodes for the letter
        surroundingNodes.(letter) = surroundingNodesForLetter;
    end
end
