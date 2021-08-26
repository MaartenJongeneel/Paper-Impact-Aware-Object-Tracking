function x3 = xprod(x1,x2)
% Group operation on SO(3)xR3xR3xR3
%
% INPUTS:    x1,x2   : states: elements of the group SO(3)xR3xR3xR3
%
% OUTPUTS:   x3      : Result of the group operation, also on
%                      SO(3)xR3xR3xR3
%% Script
R1 = x1{1};     R2 = x2{1};
o1 = x1{2};     o2 = x2{2};
v1 = x1{3};     v2 = x2{3};
w1 = x1{4};     w2 = x2{4};

%Group operation:
x3{1,1} = R1*R2;
x3{2,1} = o1+o2;
x3{3,1} = v1+v2;
x3{4,1} = w1+w2;


