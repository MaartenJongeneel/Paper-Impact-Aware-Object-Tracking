function xinv = invx(x)
% Inverse of an element of the Lie group SO(3)xR3xR3xR3
%
% INPUTS:    x       : Element of the Lie group SO(3)xR3xR3xR3
%                      represented as cell(R;o;v;w)
%
% OUTPUTS:   xinv    : Inverse of that element represented as
%                      cell(R;o;v;w)
%% Script
xinv{1,1} = inv(x{1,1});
xinv{2,1} = -x{2,1};
xinv{3,1} = -x{3,1};
xinv{4,1} = -x{4,1};