function [XNew]=CVModel(Xold,const)
% Computes the likelihood 
% INPUTS:    Xold         : cell(R;o;v;w): state of the particle 
%
% OUTPUTS:   Xnew         : cell(R;o;v;w): state of the particle 
%% Compute the state update with constant velocity model
te    = 1/60;     %Duration of one timestep (1/60 sec) [s]
    
Hold = [Xold{1} Xold{2}; zeros(1,3) 1]; %AH_B at beginning
Vold = [Xold{3};Xold{4}];               %AV_B at beginning

noise = const.Noise*randn(6,1); %Noise vector

Hnew = Hold*expm(te*hat(Vold)+sqrt(te)*hat(noise)); %Configuration update
Vnew = vee(logm((Hold)\Hnew)/te);                   %Velocity update

XNew{1} = Hnew(1:3,1:3);
XNew{2} = Hnew(1:3,4);
XNew{3} = Vnew(1:3);
XNew{4} = Vnew(4:6);
