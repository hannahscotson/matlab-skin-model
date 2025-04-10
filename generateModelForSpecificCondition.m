function generateModelForSpecificCondition(nNodes, nMass, nDoF, nTriangles, nMeshNodes, InputForce, InpM, InputSigma, v, p)
    % Medium C
    boundaryNodes = [38 39 40 41 53 58 68 75 83 92 99 115 131 147 163 172 180 187 197 202 214 215 216 217];
    InpC = 8;
    InpK = 460;
    
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
    
    % Extract the full response over time (z-displacement at all time steps for the first layer)
    model_response_over_time = x1(:, z_dof_indices); % Extract z-displacement for all time steps
    
    % Reshape each time step's displacement data to match 16x16 grid (flattened form for saving)
    num_time_steps = size(model_response_over_time, 1);
    reshaped_response = reshape(model_response_over_time', [256, num_time_steps]);  % Reshape to 256 rows (nodes), with each column representing a time step
    
    % Save the full force response over time to a .mat file
    disp(['Current directory: ', pwd]);  % This will display the current directory
    save('model_FSR_full_response.mat', 'reshaped_response', 'ts');
    
    % Print completion message
    fprintf('Completed: Medium C\n');
end
