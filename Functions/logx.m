function X = logx(x)
% Logarithmic mapping on SO(3)xR3xR3xR3
%
% INPUTS:    x      : state at SO(3)xR3xR3xR3 representend as cell(R;o;v;w)
%
% OUTPUTS:   X      : mapping to tangent space represented in R12
%
%% Script
X(1:3,1) = vee(logm(x{1}));     
X(4:6,1) = x{2};     
X(7:9,1) = x{3};     
X(10:12,1) = x{4}; 