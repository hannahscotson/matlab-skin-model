function [K,D] = BoundryNodes(K,D,nMeshNodes,nDoF,nNodes,nMass)
    % Layer 1 - Fixing nodes (Borders)
    % Bottom layer
    nMassBtmFixed = nMeshNodes(1);
    K(1:nMassBtmFixed*nDoF,:) = 0;  K(:,1:nMassBtmFixed*nDoF) = 0;
    
    % Left-side layer
    nMassLftFixed = 1:nMeshNodes(1)*nDoF:nMass*nDoF;
    nMassLftFixed = nMassLftFixed(1:nMeshNodes(2));
    nMassLftFixed = [nMassLftFixed nMassLftFixed+1  nMassLftFixed+2];
    nMassLftFixed = sort(nMassLftFixed);
    K(nMassLftFixed,:) = 0;  K(:,nMassLftFixed) = 0;
    
    % Right-side layer
    nMassRgtFixed = nMeshNodes(1)*nDoF:nMeshNodes(1)*nDoF:nMass*nDoF;
    nMassRgtFixed = nMassRgtFixed(1:nMeshNodes(2));
    nMassRgtFixed = [nMassRgtFixed nMassRgtFixed-1 nMassRgtFixed-2];
    nMassRgtFixed = sort(nMassRgtFixed);
    K(nMassRgtFixed,:) = 0;  K(:,nMassRgtFixed) = 0;
    
    % Top layer
    nNodes_FirLay = nMass - (nNodes.SecLay+nNodes.TrdLay);
    nMassTopFixed = (nNodes_FirLay-nMeshNodes(1))+1 : nNodes_FirLay;
    nMassTopFixed = nMassTopFixed * nDoF;
    nMassTopFixed = [nMassTopFixed nMassTopFixed-1 nMassTopFixed-2];
    nMassTopFixed = sort(nMassTopFixed);
    K(nMassTopFixed,:) = 0;  K(:,nMassTopFixed) = 0;
    
    clear nMassBtmFixed nMassLftFixed nMassRgtFixed nMassTopFixed
    
    % Layer 2 - Fixing nodes (Borders)
    % Bottom layer
    nMassBtmFixed = ((nNodes_FirLay+1)*nDoF)-(nDoF-1) : (nNodes_FirLay + (nMeshNodes(1)-1))*nDoF;
    K(nMassBtmFixed,:) = 0;  K(:,nMassBtmFixed) = 0;
    
    % Left-side layer
    nMassLftFixed = ((nNodes_FirLay+1)*nDoF)-(nDoF-1):(nMeshNodes(1)-1)*nDoF:(nMass-nNodes.TrdLay)*nDoF;
    nMassLftFixed = [nMassLftFixed nMassLftFixed+1  nMassLftFixed+2];
    nMassLftFixed = sort(nMassLftFixed);
    K(nMassLftFixed,:) = 0;  K(:,nMassLftFixed) = 0;
    
    % Right-side layer
    nMassRgtFixed = ((nNodes_FirLay+(nMeshNodes(1)-1))*nDoF)-(nDoF-1) : ...
        (nMeshNodes(2)-1)*nDoF:(nMass-nNodes.TrdLay)*nDoF;
    nMassRgtFixed = [nMassRgtFixed nMassRgtFixed+1  nMassRgtFixed+2];
    nMassRgtFixed = sort(nMassRgtFixed);
    K(nMassRgtFixed,:) = 0;  K(:,nMassRgtFixed) = 0;
    
    % Top layer
    nNodesTemp = (nMass-nNodes.TrdLay) - (nMeshNodes(1)-1);
    nMassTopFixed = nNodesTemp+1:(nMass-nNodes.TrdLay);
    nMassTopFixed = nMassTopFixed * nDoF;
    nMassTopFixed = [nMassTopFixed nMassTopFixed-1 nMassTopFixed-2];
    nMassTopFixed = sort(nMassTopFixed);
    K(nMassTopFixed,:) = 0;  K(:,nMassTopFixed) = 0;
    
    clear nMassBtmFixed nMassLftFixed nMassRgtFixed nMassTopFixed
    
    % Layer 3 - Fixing all nodes (Borders)
    % Bottom layer
    nMassBtmFixed = ((nNodes_FirLay+nNodes.SecLay+1)*nDoF)-(nDoF-1) : (nMass*nDoF);
    K(nMassBtmFixed,:) = 0;  K(:,nMassBtmFixed) = 0;
    D(nMassBtmFixed,:) = 0;  D(:,nMassBtmFixed) = 0;
    
    clear nMassBtmFixed

end