function h = expH(H)
% Exponential mapping onto the Lie group SO(3)xR3
%
% INPUTS:    H       : Element of the Lie algebra so(3)xR3 represented in
%                      R6
%
% OUTPUTS:   h       : Element of the Lie group SO(3)xR3 represented
%                      as cell(R;o)
%% Script
h{1,1} = expm(hat(H(1:3))); %Writes H(1:3) as element of so(3) and computes the exponential to SO(3)
h{2,1} = [H(4);H(5);H(6)];  %Origin
