function [t,p,nNodes,nT_FirLay] = MeshTriangles(nMeshNodes,nDoF,dMass,dLayer)
    
    %% First layer connections
    m=nMeshNodes(1); n=nMeshNodes(2); % includes boundary nodes, mesh spacing 1/(m-1) and 1/(n-1)
    [x,y]=ndgrid((0:m-1)/(m-1),(0:n-1)/(n-1)); % matlab forms x and y lists
    p=[x(:),y(:)]; % N by 2 matrix listing x,y coordinates of all N=mn nodes
    t=[1,2,m+2; 1,m+2,m+1]; % 3 node numbers for two triangles in first square
    t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');
    % now t lists 3 node numbers of 2(m-1) triangles in the first mesh row
    t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
    % number of nodes in first layer
    nNodes.FirLay = numel(p(:,1));
    nT_FirLay = numel(t(:,1));
    
    %% Second layer connections
    m2=nMeshNodes(1)-1; n2=nMeshNodes(2)-1; % includes boundary nodes, mesh spacing 1/(m-1) and 1/(n-1)
    t2=[1,2,m2+2; 1,m2+2,m2+1]; % 3 node numbers for two triangles in first square
    t2=kron(t2,ones(m2-1,1))+kron(ones(size(t2)),(0:m2-2)');
    % now t lists 3 node numbers of 2(m-1) triangles in the first mesh row
    t2=kron(t2,ones(n2-1,1))+kron(ones(size(t2)),(0:n2-2)'*m2);
    % Renumbering Node indexes
    t2 = t2 + nNodes.FirLay;
    % ofsetting the second layer to be in center with first layer 
    % this is done to achieve the pyramid structure
    for idx = 1:numel(t(:,1))/2
        idxTemp = [t(idx,1) t(idx,3)];
        p_SecLay(idx,1) = p(idxTemp(1),1) + ((p(idxTemp(2),1) - p(idxTemp(1),1))/2);
        p_SecLay(idx,2) = p(idxTemp(1),2) + ((p(idxTemp(2),2) - p(idxTemp(1),2))/2);
    end
    % concatenate position vector for 1st & 2nd layers of masses & triangles
    p = [p ; p_SecLay];
    t = [t;t2];
    % number of nodes in second layer
    tempY = y(1:end-1,1:end-1);
    nNodes.SecLay = numel(tempY);
    
    %% 3rd Layer Connections (Triangles in Third Layer)
    % generate mesh of T=2(m-1)(n-1) right triangles in unit square
    m3=m2-1; n3=n2-1; % includes boundary nodes, mesh spacing 1/(m-1) and 1/(n-1)
    t3=[1,2,m3+2; 1,m3+2,m3+1]; % 3 node numbers for two triangles in first square
    t3=kron(t3,ones(m3-1,1))+kron(ones(size(t3)),(0:m3-2)');
    % now t lists 3 node numbers of 2(m-1) triangles in the first mesh row
    t3=kron(t3,ones(n3-1,1))+kron(ones(size(t3)),(0:n3-2)'*m3);
    % Renumbering Node indexes
    t3 = t3 + nNodes.FirLay + nNodes.SecLay;
    % ofsetting the third layer to be in center with second layer 
    % this is done to achieve the pyramid structure
    for idx = 1:numel(t2(:,1))/2
        idxTemp = [t2(idx,1) t2(idx,3)];
        p_TrdLay(idx,1) = p(idxTemp(1),1) + ((p(idxTemp(2),1) - p(idxTemp(1),1))/2);
        p_TrdLay(idx,2) = p(idxTemp(1),2) + ((p(idxTemp(2),2) - p(idxTemp(1),2))/2);
    end
    % concatenate position vector of 1st, 2nd & 3rd layers masses & triangles
    p = [p ; p_TrdLay];
    t = [t;t3];
    % number of nodes in third layer
    nNodes.TrdLay = numel(p(:,1)) - (nNodes.FirLay + nNodes.SecLay);
    
    %% Defining the distance between the layers (z-axis)
    if nDoF == 3
        tempZ1 = ones(nNodes.FirLay,1) * (dLayer*2);
        tempZ2 = ones(nNodes.SecLay,1) * (dLayer);
        tempZ3 = ones(nNodes.TrdLay,1) * 0;
        p(:,3) = [tempZ1;tempZ2;tempZ3];
    end
    
    %% Adjusting position vector based on the distance between masses
    SameLayDist = nMeshNodes(1) * dMass;
    p(:,1) = p(:,1)* SameLayDist;
    p(:,2) = p(:,2)* SameLayDist;
    
    %% Defining the traingles between layers 
    % 1st to 2nd layer connections.
    nT_SecLay = (m*n)-n; % number of triangles in second layer
    tempT = zeros(nT_SecLay,3);
    idxNodeEachRow = 1:m:(m*n);
    if numel(idxNodeEachRow)>2
        u = repelem(idxNodeEachRow,2);
        u(1) = []; u(end) = [];
    end
    % co-ordinates on x- , y-axis
    cnt = 0;
    for j = u
        for i = 1:m-1
            cnt = cnt +1;
            tempT(cnt,1:2) = (j-1) + [i i+1];
        end
    end
    % co-ordinates on z-axis
    tempVal = nNodes.FirLay+1:(nNodes.FirLay+nNodes.SecLay);
    tempVal = reshape(tempVal,m-1,n-1);
    tempVect = [];
    for i = 1:n-1
        tempVect = [tempVect;tempVal(:,i);tempVal(:,i)];
    end
    tempT(:,3) = tempVect;
    % concatenate triangles
    t = [t ; tempT];
    clear tempVal tempVect
    
    % 2nd to 3rd layer connections.
    nT_TrdLay = (m2*n2)-n2; % number of triangles in third layer
    tempT = zeros(nT_TrdLay,3);
    idxNodeEachRow2 = nNodes.FirLay + (1:m2:(m2*n2));
    if numel(idxNodeEachRow2)>2
        u2 = repelem(idxNodeEachRow2,2);
        u2(1) = []; u2(end) = [];
    end
    % co-ordinates on x- , y-axis
    cnt = 0;
    for j = u2
        for i = 1:m2-1
            cnt = cnt +1;
            tempT(cnt,1:2) = (j-1) + [i i+1];
        end
    end
    % co-ordinates on z-axis
    tempVal = (nNodes.FirLay+nNodes.SecLay)+1:(nNodes.FirLay+nNodes.SecLay+nNodes.TrdLay);
    tempVal = reshape(tempVal,m2-1,n2-1);
    tempVect = [];
    for i = 1:n2-1
        tempVect = [tempVect;tempVal(:,i);tempVal(:,i)];
    end
    tempT(:,3) = tempVect;
    % concatenate triangles
    t = [t ; tempT];
    clear tempVal tempVect

end

