function [xPred,PxxPred,xPredSig,XPredSig,sigN,w,nsp]=ukf(xEst,PEst,Q,R,alpha,beta,kappa,fname,const)
% Compute one full step of the unscented Kalman filter
%
% INPUTS:    xEst    : state mean estimate at time k: 4x1 cell
%            PEst    : state covariance at time: k 12x12 matrix
%            Q       : process noise covariance at time k
%            R       : measurement noise covariance at k+1
%            alpha   : sigma point scaling parameter.
%            beta    : higher order error scaling parameter.
%            kappa   : scalar tuning parameter 1.
%
% OUTPUTS:   xEst    : updated estimate of state mean at time k+1
%            PEst    : updated state covariance at time k+1
%            xPred   : prediction of state mean at time k+1
%            PxxPred : prediction of state covariance at time k+1
%% Script
fname = str2func(fname);

% Calculate the dimensions of the problem 
sDim = 12; %State dimension
vDim = 12; %Process noise dimension
wDim = 6;  %Measurement noise dimension

%Augment the state mean and covariance with the noise parameters.
PA=[                PEst, zeros(sDim,vDim), zeros(sDim,wDim);
    zeros(vDim,sDim),                    Q, zeros(vDim,wDim);
    zeros(wDim,sDim), zeros(wDim,vDim),                    R];
 
xA=[zeros(sDim,1);zeros(vDim,1);zeros(wDim,1)];

%Calculate the sigma points and their weights using the Scaled Unscented
%Transformation (on tangent space of augmented state)
[sigX,sigM,sigN,w,nsp]=CompSigmaPnts(xA,PA,alpha,beta,kappa,sDim,vDim,wDim); 

for ii = 1:nsp % For each sigma point
%Map the sigma point of the state to the group
sigx = xprod(xEst,expx(sigX(:,ii))); 

%Compute the predicted sigma points on the group: f(x)*exp(n)
xPredSig(:,ii) = xprod(feval(fname,sigx,const),expx(sigM(:,ii)));   
xmean = xPredSig(:,1); %Take zeroth sigma point as initial mean on group
XPredSig(:,ii) = logx(xprod(invx(xmean),xPredSig(:,ii)));
end 

% Compute the mean in the tangent space
Xmean = XPredSig*w(1:nsp)';

% Check update and compute new mean
while norm(Xmean) > 1e-4 %If mean is not at the origin
    xmean = xprod(xmean,expx(Xmean)); %New mean on Lie group
    for ii = 1:nsp
        XPredSig(:,ii) = logx(xprod(invx(xmean),xPredSig(:,ii)));
    end
    Xmean = XPredSig*w(1:nsp)';
end

wnew=(flipud(w'))'; %To ensure zeroth weight of covariance is at index 1

%Final mean and covariance of the predicted state
PxxPred = (wnew(1:nsp).*XPredSig)*XPredSig';
xPred = xmean;        
