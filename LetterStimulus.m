function [f] = NewLetterStimulus(nMass, nDoF, InputForce, boundaryNodes)

    % Force application using Indent case
    f = zeros(nMass * nDoF, 1); % Initialize the force vector
    idxMass = 1:nDoF:nMass*nDoF; % Get starting indices for all nodes
    selectedNodes = boundaryNodes; % Get selected boundary nodes
    
    % Convert node indices to force vector indices (z-component)
    idxFz = idxMass(selectedNodes) + 2;  % Get z-axis force indices
    
    % Apply force to the correct DoFs
    f(idxFz) = InputForce;
    
    % DEBUGGING
    disp('Selected Nodes:'); disp(selectedNodes);
    disp('Force Indices (Z-axis):'); disp(idxFz);
    disp('Applied Forces:'); disp(f(f ~= 0));
end