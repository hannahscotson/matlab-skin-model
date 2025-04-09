function generateModelForSpecificCondition(nNodes, nMass, nDoF, nTriangles, nMeshNodes, InputForce, InpM, InputSigma, v, p, targetLetter, targetHardness)
% Letters and hardness levels
letters = {'C', 'D', 'Q'};
hardness_levels = {'Soft', 'Medium', 'Hard'};

% Find the indices for the target letter and hardness
letter_idx = find(strcmp(letters, targetLetter));
hardness_idx = find(strcmp(hardness_levels, targetHardness));

% Check if the target combination exists
if isempty(letter_idx) || isempty(hardness_idx)
    error('Invalid letter or hardness level.');
end

% Define parameters for the chosen letter
switch upper(targetLetter)
    case 'C'
        boundaryNodes = [38 39 40 41 53 58 68 75 83 92 99 115 131 147 163 172 180 187 197 202 214 215 216 217];
    case 'D'
        boundaryNodes = [36 37 38 39 40 41 52 58 68 75 84 92 100 108 116 124 132 140 148 156 164 172 180 187 196 202 212 213 214 215 216 217];
    case 'Q'
        boundaryNodes = [38 39 40 41 42 46 53 59 61 68 76 83 91 93 99 109 115 125 131 141 147 157 163 173 180 188 197 203 214 215 216 217 218];
end

% Adjust parameters based on hardness level
switch targetHardness
    case 'Soft'
        InpC = 11;
        InpK = 160;
    case 'Medium'
        InpC = 8;
        InpK = 460;
    case 'Hard'
        InpC = 7;
        InpK = 690;
end

% Defining the values of spring/damper for individual layers
Kval.L1 = InpK;  Kval.L2 = InpK; Kval.L3 = InpK;
Cval.L1 = InpC;  Cval.L2 = InpC; Cval.L3 = InpC;
[valK,valC] = SpringVals(v,Kval,Cval,nNodes,nMass);

% Mass matrix
valM = ones(nMass*nDoF,1)* InpM;
M = diag(valM);

% Compute stiffness & damping matrices
[K,D] = stiffnessMatrix(nMass,nTriangles,valK,valC,v,p,nDoF);

% Defining boundary nodes
[K,D] = BoundryNodes(K,D,nMeshNodes,nDoF,nNodes,nMass);

% Call the LetterStimulus function to generate force for the selected letter
[f] = LetterStimulus(nMass, nDoF, InputForce, boundaryNodes);

% Simulation time
ts = 0.18:1e-5:0.245;

% Gaussian force profile in time
[Fx,FxDyy,tFx] = GaussF(nMass,nDoF,ts,InputSigma,f);

% Initial position vector
u = zeros(nMass*nDoF,1);
u(1) = 1e-3;  % Small initial displacement (e.g., 0.001 m or 1 mm)

% Initial velocity vector
udot = zeros(nMass*nDoF,1);

% Initial vector
uo = [u ; udot];

% Solve Mass-Spring-Damper system
[t, x] = ode45(@(ts, uo) MassSpringDamper(ts, uo, K, D, M, nMass, FxDyy, tFx, nDoF), ts, uo);

% Position & velocity
x1 = x(:,1:nMass*nDoF); % Extract position data
ts = t; % Extract time steps

% Extracting the z-displacement (3rd DoF of each node) for the first layer nodes
first_layer_nodes = 1:256; % Indices of first layer nodes
z_dof_indices = (first_layer_nodes - 1) * 3 + 3; % Z-axis DoF indices

% Find the timestep where the max absolute displacement occurs (across all nodes)
[~, max_timestep] = max(max(abs(x1(:, z_dof_indices)), [], 2)); % Find max frame

% Extract data at this timestep and reshape to 16x16
model_response_max = reshape(x1(max_timestep, z_dof_indices), [16, 16]);

% Store the response in FSR format (flattened 16x16 → 1x256)
new_trial_data = reshape(model_response_max, [1, 256]);  % Store current trial response

% Save the unnormalised data to a file
save('model_FSR_unnormalised.mat', 'new_trial_data');

% Print completion message
fprintf('Completed: %s - %s\n', targetLetter, targetHardness);

end
