function v=MassSpringDamper(ts,x,K,D,M,nMass,FxDyy,tFx,nDoF)
    % DEBUGGING
        disp(['Solving at t = ', num2str(ts)]);
        pause(0.1); % Small delay to see output
        % Your function code...
    
    A = [zeros(nMass*nDoF) eye(nMass*nDoF);-inv(M)*(K) -inv(M)*D];
    
    b = [zeros(nMass*nDoF);inv(M)];
    
    Fx = interp1(tFx,FxDyy',ts); Fx = Fx';
    
    v = (A*x) + (b*Fx);

end

% ts: The current time at which the state is being evaluated.
% x: The state vector containing positions and velocities of all mass points.
% K: The stiffness matrix of the system.
% D: The damping matrix of the system.
% M: The mass matrix of the system.
% nMass: The number of mass points in the system.
% FxDyy: A matrix (or vector) representing the external forces applied to the system.
% tFx: The time vector corresponding to the force data.
% nDoF: The number of degrees of freedom per mass.

% v: The time derivative of the state vector, which represents the rates of
% change of positions and velocities.

% matrix A:
% The top-left block is a zero matrix of size nMass×nDoF, which corresponds to
%     the derivatives of the positions (i.e., the velocities).
% The top-right block is the identity matrix of size nMass×nDoF, which means
%     the derivatives of the positions with respect to time yield the current velocities.
% The bottom-left block is the result of -M^(-1)K, which relates to the accelerations
%     derived from Hooke-s Law (stiffness).
% The bottom-right block is -M^(-1)D, which relates to the damping forces acting on the system.

% matrix b:
% The top part is a zero vector of size nMass×nDoF (for position derivatives).
% The bottom part is the inverse of the mass matrix M (for the influence of forces on acceleration).

