function p = LieProb(state,mean,cov)
% Compute probability: normal distribution on Lie group
% INPUTS:    state        : cell(R;o;v;w): state of the particle 
%            mean         : cell(R;o;v;w): mean of the distribution
%            cov          : 12x12 matrix : covariance of the distribution
%
% OUTPUTS:   p            : Probability 
%% Compute the probability
Beta = 1/(sqrt((2*pi)^12*det(cov)));    %Scaling parameter

zeta = logx(xprod(invx(mean),state));   %Something without name ?

p = Beta*exp(-0.5*zeta'*inv(cov)*zeta); %Normal distribution 