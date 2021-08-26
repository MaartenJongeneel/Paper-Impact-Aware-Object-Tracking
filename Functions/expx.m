function x = expx(X)
% Exponential mapping onto the Lie group SO(3)xR3xR3xR3
% 
% INPUTS:    X       : Element of the Lie algebra so(3)xR3xR3xR3
%                      represented in R12
% 
% OUTPUTS:   x       : Element of the Lie group SO(3)xR3xR3xR3 represented
%                      as cell(R;o;v;w)
%% Script
x{1,1} = expm(hat(X(1:3)));   %Writes X(1:3) as element of so(3) and compute the exponential to SO(3)
x{2,1} = [X(4);X(5);X(6)];    %Origin
x{3,1} = [X(7);X(8);X(9)];    %Linear velocity
x{4,1} = [X(10);X(11);X(12)]; %Angular velocity