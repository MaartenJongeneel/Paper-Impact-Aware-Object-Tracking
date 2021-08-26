function [sigX,sigM,sigN,wPts,nPts] = CompSigmaPnts(x,P,alpha,beta,kappa,sDim,vDim,wDim)
% Compute one full step of the unscented Kalman filter
%
% INPUTS:    x       : state mean estimate at time k
%            P       : state covariance at time k
%            alpha   : sigma point scaling parameter.
%            beta    : higher order error scaling parameter.
%            kappa   : scalar tuning parameter 1.
%            sDim    : State dimension
%            vDim    : Process noise dimension
%            wDim    : Measurement noise dimension
%
% OUTPUTS:   sigX    : the sigma points of the state
%            sigM    : the sigma points of the process noise
%            sigN    : the sigam points of the measurement noise
%            wPts    : the weights on the points
%            nPts    : the number of points
%% Script
%Number of sigma points and scaling terms
n    = size(x(:),1); %Size of the augmented state
nPts = 2*n+1;        %Using symmetric SUT: 2 points for each dim + mean

%Compute Lambda according to scaling parameters
lambda = alpha^2*(n+kappa)-n;

%Compute matrix square root of weighted covariance matrix (Cholesky decom.)
Psqrtm=(sqrtm((n+lambda)*P));  

%Array of the sigma points
xPts=[zeros(size(P,1),1) Psqrtm -Psqrtm];

%Add mean back in
xPts = xPts + repmat(x,1,nPts);  

%Array of the weights for each sigma point
wPts=[lambda 0.5*ones(1,nPts-1) 0]/(n+lambda);

%Now calculate the zero'th covariance term weight
wPts(nPts+1) = wPts(1) + (1-alpha^2) + beta;

sigX = xPts(1:sDim,:);                     %State sigma points
sigM = xPts(sDim+1:sDim+vDim,:);           %Process noise sigma points
sigN = xPts(sDim+vDim+1:sDim+vDim+wDim,:); %Measurement noise sigma points
