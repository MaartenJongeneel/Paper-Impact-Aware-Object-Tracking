function hinv = invH(h)
% Inverse of an element of the Lie group SO(3)xR3
%
% INPUTS:    h       : Element of the Lie group SO(3)xR3
%                      represented as cell(R;o)
%
% OUTPUTS:   hinv    : Inverse of that element represented as
%                      cell(R;o;)
%% Script
hinv{1,1} = inv(h{1,1});
hinv{2,1} = -h{2,1};