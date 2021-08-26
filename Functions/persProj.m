function uv=persProj(xyz,K)
% This function computes the 2D points in the image plane from the 3D world
% coordinates xyz. (Perspective Projection)
%
% INPUTS:    xyz       : [x y z] or [x; y; z] coordinates in the 3D world
%            K         : Intrinsic camera parameter matrix
%
% OUTPUTS:   uv        : [u;v] coordinates of the points in the image plane
%% Compute the points
uv=K*(xyz./xyz(3,:));  % uv_point=K*(xyz_point/z);
