function [K,D] = stiffnessMatrix(nMass,nTriangles,valK,valC,v,p,nDoF)
% Computation of stiffness and damping matrices

% index of each all masses in the skin
idxMass = 1:nDoF:nMass*nDoF;
% Variable definition
K = zeros(nMass*nDoF);
D = zeros(nMass*nDoF);

% Computing the stiffness and damping matrix for each triangle in the skin
for T = 1:nTriangles

    tK1 = zeros(nMass*nDoF); tK2 = zeros(nMass*nDoF); tK3 = zeros(nMass*nDoF);
    tD1 = zeros(nMass*nDoF); tD2 = zeros(nMass*nDoF); tD3 = zeros(nMass*nDoF);

    % Mass for given triangle
    Mi = v(T,1); Mj = v(T,2); Mk = v(T,3);

    % Phi being computed based on distance bewtween masses
    % Mi -> Mj
    xi = p(Mi,1); yi = p(Mi,2); zi = p(Mi,3); xj = p(Mj,1); yj = p(Mj,2); zj = p(Mj,3);
    L = sqrt((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2);
    c1 = (xj-xi)/L; s1 = (yj-yi)/L; t1 = (zj-zi)/L;
    clear L xi yi zi zj xj yj

    % Mj -> Mk
    xi = p(Mj,1); yi = p(Mj,2); zi = p(Mj,3); xj = p(Mk,1); yj = p(Mk,2); zj = p(Mk,3);
    L = sqrt((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2);
    c2 = (xj-xi)/L; s2 = (yj-yi)/L; t2 = (zj-zi)/L;
    clear L xi yi zi zj xj yj

    % Mi -> Mk
    xi = p(Mi,1); yi = p(Mi,2); zi = p(Mi,3); xj = p(Mk,1); yj = p(Mk,2); zj = p(Mk,3);
    L = sqrt((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2);
    c3 = (xj-xi)/L; s3 = (yj-yi)/L; t3 = (zj-zi)/L;
    clear L xi yi zi zj xj yj


    % Spring stiffness values for the springs in the given triangle
    index = and(any(valK(:,1:2) == Mi,2), any(valK(:,1:2) == Mj,2));
    k1 = valK(index == 1,3);
    clear index
    index = and(any(valK(:,1:2) == Mj,2), any(valK(:,1:2) == Mk,2));
    k2 = valK(index == 1,3);
    clear index
    index = and(any(valK(:,1:2) == Mk,2), any(valK(:,1:2) == Mi,2));
    k3 = valK(index == 1,3);
    clear index
    % damper coefficient for the dampers in the given triangle
    index = and(any(valC(:,1:2) == Mi,2), any(valC(:,1:2) == Mj,2));
    d1 = valC(index == 1,3);
    clear index
    index = and(any(valC(:,1:2) == Mj,2), any(valC(:,1:2) == Mk,2));
    d2 = valC(index == 1,3);
    clear index
    index = and(any(valC(:,1:2) == Mk,2), any(valC(:,1:2) == Mi,2));
    d3 = valC(index == 1,3);
    clear index

    %% stiffness matrix
    % Mi -> Mj
    A1 = [k1 -k1; -k1 k1];
    C1 = [d1 -d1; -d1 d1];
    B1 = [c1^2 c1*s1 c1*t1; c1*s1 s1^2 s1*t1; c1*t1 s1*t1 t1^2];
    K1 = kron(A1,B1);
    D1 = kron(C1,B1);

    tK1(idxMass(Mi):idxMass(Mi)+2,idxMass(Mi):idxMass(Mi)+2) = K1(1:3,1:3);
    tK1(idxMass(Mj):idxMass(Mj)+2,idxMass(Mj):idxMass(Mj)+2) = K1(4:6,4:6);
    tK1(idxMass(Mi):idxMass(Mi)+2,idxMass(Mj):idxMass(Mj)+2) = K1(1:3,4:6);
    tK1(idxMass(Mj):idxMass(Mj)+2,idxMass(Mi):idxMass(Mi)+2) = K1(4:6,1:3);

    tD1(idxMass(Mi):idxMass(Mi)+2,idxMass(Mi):idxMass(Mi)+2) = D1(1:3,1:3);
    tD1(idxMass(Mj):idxMass(Mj)+2,idxMass(Mj):idxMass(Mj)+2) = D1(4:6,4:6);
    tD1(idxMass(Mi):idxMass(Mi)+2,idxMass(Mj):idxMass(Mj)+2) = D1(1:3,4:6);
    tD1(idxMass(Mj):idxMass(Mj)+2,idxMass(Mi):idxMass(Mi)+2) = D1(4:6,1:3);


    % Mj -> Mk
    A2 = [k2 -k2; -k2 k2];
    C2 = [d2 -d2; -d2 d2];
    B2 = [c2^2 c2*s2 c2*t2; c2*s2 s2^2 s2*t2; c2*t2 s2*t2 t2^2];
    K2 = kron(A2,B2);
    D2 = kron(C2,B2);

    tK2(idxMass(Mj):idxMass(Mj)+2,idxMass(Mj):idxMass(Mj)+2) = K2(1:3,1:3);
    tK2(idxMass(Mk):idxMass(Mk)+2,idxMass(Mk):idxMass(Mk)+2) = K2(4:6,4:6);
    tK2(idxMass(Mj):idxMass(Mj)+2,idxMass(Mk):idxMass(Mk)+2) = K2(1:3,4:6);
    tK2(idxMass(Mk):idxMass(Mk)+2,idxMass(Mj):idxMass(Mj)+2) = K2(4:6,1:3);

    tD2(idxMass(Mj):idxMass(Mj)+2,idxMass(Mj):idxMass(Mj)+2) = D2(1:3,1:3);
    tD2(idxMass(Mk):idxMass(Mk)+2,idxMass(Mk):idxMass(Mk)+2) = D2(4:6,4:6);
    tD2(idxMass(Mj):idxMass(Mj)+2,idxMass(Mk):idxMass(Mk)+2) = D2(1:3,4:6);
    tD2(idxMass(Mk):idxMass(Mk)+2,idxMass(Mj):idxMass(Mj)+2) = D2(4:6,1:3);

    % Mk -> Mi
    A3 = [k3 -k3; -k3 k3];
    C3 = [d3 -d3; -d3 d3];
    B3 = [c3^2 c3*s3 c3*t3; c3*s3 s3^2 s3*t3; c3*t3 s3*t3 t3^2];
    K3 = kron(A3,B3);
    D3 = kron(C3,B3);

    tK3(idxMass(Mk):idxMass(Mk)+2,idxMass(Mk):idxMass(Mk)+2) = K3(1:3,1:3);
    tK3(idxMass(Mi):idxMass(Mi)+2,idxMass(Mi):idxMass(Mi)+2) = K3(4:6,4:6);
    tK3(idxMass(Mk):idxMass(Mk)+2,idxMass(Mi):idxMass(Mi)+2) = K3(1:3,4:6);
    tK3(idxMass(Mi):idxMass(Mi)+2,idxMass(Mk):idxMass(Mk)+2) = K3(4:6,1:3);

    tD3(idxMass(Mk):idxMass(Mk)+2,idxMass(Mk):idxMass(Mk)+2) = D3(1:3,1:3);
    tD3(idxMass(Mi):idxMass(Mi)+2,idxMass(Mi):idxMass(Mi)+2) = D3(4:6,4:6);
    tD3(idxMass(Mk):idxMass(Mk)+2,idxMass(Mi):idxMass(Mi)+2) = D3(1:3,4:6);
    tD3(idxMass(Mi):idxMass(Mi)+2,idxMass(Mk):idxMass(Mk)+2) = D3(4:6,1:3);

    K = K + tK1; K = K + tK2; K = K + tK3;
    D = D + tD1; D = D + tD2; D = D + tD3;

    clear tK1 tK2 tK3 s1 s2 s3 c1 c2 c3 t1 t2 t3 k1 k2 k3 Mi Mj Mk A1 B1 K1 A2 B2 K2 A3 B3 K3
    clear tD1 tD2 tD3 d1 d2 d3 C1 D1 C2 D2 C3 D3

end
end

