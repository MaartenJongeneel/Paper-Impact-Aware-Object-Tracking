function H3 = Hprod(H1,H2)
% Group operation on SO(3)xR3
%
% INPUTS:    H1,H2   : Poses, elements of the group SO(3)xR3
%
% OUTPUTS:   H3      : Result of the group operation, also on
%                      SO(3)xR3
%% Script
R1 = H1{1};     R2 = H2{1};
o1 = H1{2};     o2 = H2{2};

%Group operation:
H3{1,1} = R1*R2;
H3{2,1} = o1+o2;
