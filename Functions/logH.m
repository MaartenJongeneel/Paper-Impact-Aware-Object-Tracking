function H = logH(h)
% Logarithmic mapping on SO(3)xR3
%
% INPUTS:    h      : state at SO(3)xR3 representend as cell(R;o)
%
% OUTPUTS:   H      : mapping to tangent space represented in R6
%
%% Script
H(1:3,1) = vee(logm(h{1}));     
H(4:6,1) = h{2};     
 