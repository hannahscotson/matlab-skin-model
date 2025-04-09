function [Fx,FxDyy,tFx] = GaussF(nMass,nDoF,ts,InputSigma,f)
% gaussian force profile
mu  = 0.2; tFx = ts; fGain = 0.5e5;
Fx = zeros(nMass*nDoF,numel(ts));
FxDyy = zeros(nMass*nDoF,numel(ts));
for i = 1:nMass*nDoF
    Fx(i,:)     = f(i) * fGain * (gaussmf(tFx,[InputSigma,mu]));
    FxDyy(i,:)  = del2(Fx(i,:));
end

% Debugging: Check for NaN values
    if any(isnan(Fx(:))) || any(isnan(FxDyy(:)))
        error('NaN detected in force matrix Fx or FxDyy!');
    end
    
end