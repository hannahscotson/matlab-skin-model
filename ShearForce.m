function Fs = ShearForce(x1,nSprings,valK,nMass,t,nDoF,p)

% Fs -> Shear Force
% K  -> Stiffness
% x  -> Rate of change from equilibrium position

idxX = 1:nDoF:nMass*nDoF;
idxY = 2:nDoF:nMass*nDoF;
idxZ = 2:nDoF:nMass*nDoF;

Fs = zeros(numel(t),nSprings);

for i = 1:nSprings
    
    M1 = valK(i,1);
    M2 = valK(i,2);
    
    % Equilibrium
    EqM1x = p(M1,1); EqM1y = p(M1,2); EqM1z = p(M1,3);
    EqM2x = p(M2,1); EqM2y = p(M2,2); EqM2z = p(M2,3);
    
    for n = 1:numel(t)
        M1PosX = x1(n,idxX(M1));    eqM1PosX = EqM1x;
        M1PosY = x1(n,idxY(M1));    eqM1PosY = EqM1y;
        M1PosZ = x1(n,idxZ(M1));    eqM1PosZ = EqM1z;
        FnM1(n,i) = (sqrt((eqM1PosX - M1PosX)^2 + (eqM1PosY - M1PosY)^2 + (eqM1PosZ - M1PosZ)^2));
        
        M2PosX = x1(n,idxX(M2));    eqM2PosX = EqM2x;
        M2PosY = x1(n,idxY(M2));    eqM2PosY = EqM2y;
        M2PosZ = x1(n,idxZ(M2));    eqM2PosZ = EqM2z;
        FnM2(n,i) = (sqrt((eqM2PosX - M2PosX)^2 + (eqM2PosY - M2PosY)^2 + (eqM2PosZ - M2PosZ)^2));
    end
    
    Fs(:,i) = FnM1(:,i) - FnM2(:,i);

end

for i = 1:nSprings; Fs(:,i) = Fs(:,i) - Fs(1,i); end

end