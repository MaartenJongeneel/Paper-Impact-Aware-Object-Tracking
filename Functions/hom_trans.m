function points=hom_trans(pts,state)
% This function translates/rotates the unit points to the curent state of
% the cuboid.
%
% INPUTS:    pts       : [x; y; z] coordinates in the 3D world (unit points)
%            state     :  cell(R;o;v;w): state of the particle 
%
% OUTPUTS:   points    : [x; y; z] coordinates of translated points
%% Rotate/translate the point
% Rotation matrix is R = exp(w*theta) with w the first three indicies of
% the state and theta the fourth index
R = state{1};
xyz = state{2};

%Rotate the points according to rotation matrix R
points= R*pts;
%Translate the points from origin to XYZ position
points= points + repmat(xyz, 1, size(points,2));
