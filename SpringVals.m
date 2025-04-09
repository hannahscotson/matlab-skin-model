function [valK,valC] = SpringVals(v,StatKval,StatDval,nNodes,nMass)

% Springs in all triangles
for i = 1:numel(v(:,1))
    tempVal{i} = nchoosek(v(i,:),2);
end
% Concatination
tempVal = vertcat(tempVal{:});

% Unique springs (as all triangles interconnected)
[valK,~,~] = unique(tempVal,'rows');

% Removing spring repeated indexes (e.g., [a,b] = [b,a])
idxRem = [];
for i = 1:numel(valK(:,1))
    index = and(any(valK(i,1) == valK(:,2),2), any(valK(i,2) == valK(:,1),2));
    tempIdx = find(index);
    idxRem = [idxRem tempIdx];
end
valK(idxRem(1:2:end),:) = [];

% Identiying springs in different layers
tempIdx_L1 = find(valK(:,2) <= (nMass-(nNodes.SecLay+nNodes.TrdLay)));

tempIdx_L2 = find(and(valK(:,2) > (nMass-(nNodes.SecLay+nNodes.TrdLay)),...
    valK(:,2) <= (nMass-nNodes.TrdLay)));

tempIdx_L3 = find(valK(:,2) > (nMass-nNodes.TrdLay));

% Assigning the stiffness value based on layer
valK(tempIdx_L1,3) = StatKval.L1;
valK(tempIdx_L2,3) = StatKval.L2;
valK(tempIdx_L3,3) = StatKval.L3;

% Assigning the damping coefficient based on layer
valC(:,1:2) = valK(:,1:2);
valC(tempIdx_L1,3) = StatDval.L1;
valC(tempIdx_L2,3) = StatDval.L2;
valC(tempIdx_L3,3) = StatDval.L3;

end

