function plotAnimation(p,nNodes,nMass,nDoF,nMeshNodes,x1,ts)

Zoom = 0;
figure; % Open figure in a separate window
% time (rounding for printing on plot)
tp = ts-ts(1);
tp = tp *1e3; tp = round(tp); tp = tp*1e-3;

%% Initial Co-ordinates
SkinSize = nMeshNodes;
% Layer-1
Init.L1X = p(1:nNodes.FirLay,1);
Init.L1Y = p(1:nNodes.FirLay,2);
Init.L1Z = p(1:nNodes.FirLay,3);
Init.L1X = reshape(Init.L1X,SkinSize(1),SkinSize(2)); Init.L1X = Init.L1X';
Init.L1Y = reshape(Init.L1Y,SkinSize(1),SkinSize(2)); Init.L1Y = Init.L1Y';
Init.L1Z = reshape(Init.L1Z,SkinSize(1),SkinSize(2)); Init.L1Z = Init.L1Z';
% Layer-2
Init.L2X = p(nNodes.FirLay+1 : nNodes.FirLay+nNodes.SecLay,1);
Init.L2Y = p(nNodes.FirLay+1 : nNodes.FirLay+nNodes.SecLay,2);
Init.L2Z = p(nNodes.FirLay+1 : nNodes.FirLay+nNodes.SecLay,3);
Init.L2X = reshape(Init.L2X,SkinSize(1)-1,SkinSize(2)-1);
Init.L2Y = reshape(Init.L2Y,SkinSize(1)-1,SkinSize(2)-1);
Init.L2Z = reshape(Init.L2Z,SkinSize(1)-1,SkinSize(2)-1);
% Layer-3
Init.L3X = p(nNodes.FirLay+nNodes.SecLay+1 : nNodes.FirLay+nNodes.SecLay+nNodes.TrdLay,1);
Init.L3Y = p(nNodes.FirLay+nNodes.SecLay+1 : nNodes.FirLay+nNodes.SecLay+nNodes.TrdLay,2);
Init.L3Z = p(nNodes.FirLay+nNodes.SecLay+1 : nNodes.FirLay+nNodes.SecLay+nNodes.TrdLay,3);
Init.L3X = reshape(Init.L3X,SkinSize(1)-2,SkinSize(2)-2); Init.L3X = Init.L3X';
Init.L3Y = reshape(Init.L3Y,SkinSize(1)-2,SkinSize(2)-2); Init.L3Y = Init.L3Y';
Init.L3Z = reshape(Init.L3Z,SkinSize(1)-2,SkinSize(2)-2); Init.L3Z = Init.L3Z';

%% Index based on axis for each mass
idxX = 1:nDoF:nMass*nDoF;
idxY = 2:nDoF:nMass*nDoF;
idxZ = 3:nDoF:nMass*nDoF;

%% Positing values depending on mass
% Layer-1
L1_valX = zeros(nNodes.FirLay,numel(ts));
L1_valY = zeros(nNodes.FirLay,numel(ts));
L1_valZ = zeros(nNodes.FirLay,numel(ts));
for i = 1:nNodes.FirLay
    L1_valX(i,:) = x1(:,idxX(i));
    L1_valY(i,:) = x1(:,idxY(i));
    L1_valZ(i,:) = x1(:,idxZ(i));
end

% Layer-2
L2_valX = zeros(nNodes.SecLay,numel(ts));
L2_valY = zeros(nNodes.SecLay,numel(ts));
L2_valZ = zeros(nNodes.SecLay,numel(ts));
cnt = 0;
for i = nNodes.FirLay+1 : nNodes.FirLay+nNodes.SecLay
    cnt = cnt + 1;
    L2_valX(cnt,:) = x1(:,idxX(i));
    L2_valY(cnt,:) = x1(:,idxY(i));
    L2_valZ(cnt,:) = x1(:,idxZ(i));
end

% Layer-3
L3_valX = zeros(nNodes.TrdLay,numel(ts));
L3_valY = zeros(nNodes.TrdLay,numel(ts));
L3_valZ = zeros(nNodes.TrdLay,numel(ts));
cnt = 0;
for i = nNodes.FirLay+nNodes.SecLay+1 : nNodes.FirLay+nNodes.SecLay+nNodes.TrdLay
    cnt = cnt + 1;
    L3_valX(cnt,:) = x1(:,idxX(i));
    L3_valY(cnt,:) = x1(:,idxY(i));
    L3_valZ(cnt,:) = x1(:,idxZ(i));
end

%% Distance computation for defining color gradient
% Distance matric
L1_valDist = zeros(size(L1_valZ));
for i = 1:nNodes.FirLay
    for j = 2:numel(ts)
        L1_valDist(i,j) = sqrt((L1_valX(i,j) - L1_valX(i,j-1))^2 + ...
            (L1_valY(i,j) - L1_valY(i,j-1))^2 + ...
            (L1_valZ(i,j) - L1_valZ(i,j-1))^2);
    end
end
ClrMin = min(min(L1_valDist));
ClrMax = max(max(L1_valDist));

L2_valDist = zeros(size(L2_valZ));
for i = 1:nNodes.SecLay
    for j = 2:numel(ts)
        L2_valDist(i,j+1) = sqrt((L2_valX(i,j) - L2_valX(i,j-1))^2 + ...
            (L2_valY(i,j) - L2_valY(i,j-1))^2 + ...
            (L2_valZ(i,j) - L2_valZ(i,j-1))^2);
    end
end

L3_valDist = zeros(size(L3_valZ));
for i = 1:nNodes.TrdLay
    for j = 2:numel(ts)
        L3_valDist(i,j+1) = sqrt((L3_valX(i,j) - L3_valX(i,j-1))^2 + ...
            (L3_valY(i,j) - L3_valY(i,j-1))^2 + ...
            (L3_valZ(i,j) - L3_valZ(i,j-1))^2);
    end
end

% Arranging distance matrix for plotting
L1_dist = cell(numel(ts),1);
L2_dist = cell(numel(ts),1);
L3_dist = cell(numel(ts),1);

for i = 1:numel(ts)
    tL1_Z = L1_valDist(:,i);
    tL1_Z = reshape(tL1_Z,SkinSize(1),SkinSize(2));
    tL1_Z = tL1_Z';
    L1_dist{i} = tL1_Z;

    tL2_Z = L2_valDist(:,i);
    tL2_Z = reshape(tL2_Z,SkinSize(1)-1,SkinSize(2)-1);
    tL2_Z = tL2_Z';
    L2_dist{i} = tL2_Z;

    tL3_Z = L3_valDist(:,i);
    tL3_Z = reshape(tL3_Z,SkinSize(1)-2,SkinSize(2)-2);
    tL3_Z = tL3_Z';
    L3_dist{i} = tL3_Z;
end


%% Arranging mass in matrix format for plot
% Layer-1
L1_Xt = cell(numel(ts),1);
L1_Yt = cell(numel(ts),1);
L1_Zt = cell(numel(ts),1);
for i = 1:numel(ts)
    tVal = L1_valX(:,i);
    tVal = reshape(tVal,SkinSize(1),SkinSize(2)); tVal = tVal';
    L1_Xt{i} = tVal + Init.L1X;
    clear tVal

    tVal = L1_valY(:,i);
    tVal = reshape(tVal,SkinSize(1),SkinSize(2)); tVal = tVal';
    L1_Yt{i} = tVal + Init.L1Y;
    clear tVal

    tVal = L1_valZ(:,i);
    tVal = reshape(tVal,SkinSize(1),SkinSize(2)); tVal = tVal';
    L1_Zt{i} = tVal + Init.L1Z;
    clear tVal
end

% Layer-2
L2_Xt = cell(numel(ts),1);
L2_Yt = cell(numel(ts),1);
L2_Zt = cell(numel(ts),1);
for i = 1:numel(ts)
    tVal = L2_valX(:,i);
    tVal = reshape(tVal,SkinSize(1)-1,SkinSize(2)-1); tVal = tVal';
    L2_Xt{i} = tVal + Init.L2X;
    clear tVal

    tVal = L2_valY(:,i);
    tVal = reshape(tVal,SkinSize(1)-1,SkinSize(2)-1); tVal = tVal';
    L2_Yt{i} = tVal + Init.L2Y;
    clear tVal

    tVal = L2_valZ(:,i);
    tVal = reshape(tVal,SkinSize(1)-1,SkinSize(2)-1); tVal = tVal';
    L2_Zt{i} = tVal + Init.L2Z;
    clear tVal
end

% Layer-3
L3_Xt = cell(numel(ts),1);
L3_Yt = cell(numel(ts),1);
L3_Zt = cell(numel(ts),1);
for i = 1:numel(ts)
    tVal = L3_valX(:,i);
    tVal = reshape(tVal,SkinSize(1)-2,SkinSize(2)-2); tVal = tVal';
    L3_Xt{i} = tVal + Init.L3X;
    clear tVal

    tVal = L3_valY(:,i);
    tVal = reshape(tVal,SkinSize(1)-2,SkinSize(2)-2); tVal = tVal';
    L3_Yt{i} = tVal + Init.L3Y;
    clear tVal

    tVal = L3_valZ(:,i);
    tVal = reshape(tVal,SkinSize(1)-2,SkinSize(2)-2); tVal = tVal';
    L3_Zt{i} = tVal + Init.L3Z;
    clear tVal
end

%% Plot

colormap((parula));

if Zoom == 0
    for i = 1:10:numel(ts)
        surf(L1_Xt{i},L1_Yt{i},L1_Zt{i},L1_dist{i},'LineStyle','none','FaceAlpha',0.5); hold on;
        surf(L2_Xt{i},L2_Yt{i},L2_Zt{i},L2_dist{i},'LineStyle','none','FaceAlpha',0.5);
        surf(L3_Xt{i},L3_Yt{i},L3_Zt{i},L3_dist{i},'LineStyle','none','FaceAlpha',0.5);
        clim([ClrMin ClrMax/5]); hold off; % Inp5 = ClrMax / 3. Rest ClrMax/20
        title(['time ',num2str(tp(i)),' s']);
        xlim([0 5.7e-2]); ylim([0 5.7e-2]); zlim([0 1.4e-2]);
        pause(0.0001);
    end

end


if Zoom == 1
    for i = 1:10:numel(ts)
        subplot(211);
        surf(L1_Xt{i},L1_Yt{i},L1_Zt{i},L1_dist{i},'LineStyle','none','FaceAlpha',0.5);
        clim([ClrMin ClrMax/5]); % Inp5 = ClrMax / 3. Rest ClrMax/20
        title(['time ',num2str(tp(i)),' s']);
        xlim([0 15e-3]); ylim([0 15e-3]); zlim([1.75e-3 2.25e-3]);
        subplot(212);
        surf(L2_Xt{i},L2_Yt{i},L2_Zt{i},L2_dist{i},'LineStyle','none','FaceAlpha',0.5);
        clim([ClrMin ClrMax/5]); % Inp5 = ClrMax / 3. Rest ClrMax/20
        xlim([0 15e-3]); ylim([0 15e-3]); zlim([0.75e-3 1.25e-3]);
        pause(0.0001);
    end
end

end