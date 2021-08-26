function [xEst,PEst] = updateX(xPred,PxxPred,xPredSigma,XPredSigma,sigmaN,weights,nsp,Z)
for jj = 1:length(xPredSigma(1,:)) %For each sigma point
    hpred(:,jj) = Hprod(xPredSigma(1:2,jj),expH(sigmaN(:,jj))); %Hpred_j = Hpred_j*exp(Hsigma): The predicted output is perturbed with gaussian noise
    hmean = hpred(:,1);
    HPredSigma(:,jj) = logH(Hprod(invH(hmean),hpred(:,jj)));
end

%% Compute the mean
Hmean = HPredSigma*weights(1:nsp)';

% Check update and compute new mean
while norm(Hmean) > 1e-4 %Might be too computational heavy
    hmean = Hprod(hmean,expH(Hmean)); %New mean on Lie group
    for jj = 1:nsp
        HPredSigma(:,jj) = logH(Hprod(invH(hmean),hpred(:,jj)));
    end
    Hmean = HPredSigma*weights(1:nsp)';
end

wnew=flipud(weights); %To ensure zeroth weight of covariance is at index 1

%Final mean and covariance of the predicted state
PyyPred = (wnew(1:nsp).*HPredSigma)*HPredSigma';
PxyPred = (wnew(1:nsp).*XPredSigma)*HPredSigma';

inn = logH(Hprod(invH(hmean),Z)); %innovation = log(hmean^-1*Z)
xEst = xprod(xPred,expx(PxyPred/(PyyPred)*inn)); %Estimated state
PEst = PxxPred - PxyPred*(PxyPred/(PyyPred))'; %Estimated covariance