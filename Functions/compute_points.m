function [pin,pface,pout,see_face] = compute_points(box,state,K)
% This function computes the 2D points in the image plane that will be used
% to compute the inner and outer histogram of a particle
% 
% INPUTS:    box.face_centers : xyz coordinates of each face center 
%            box.Opoints   : cell array with points for each edge
%            box.Epoints   : Cell with face points for H_in for each face
%            box.Spoints  : Cell with face points for Hf for each face
%            state     : Cell(R;o;v;w): state of the particle 
%            K         : Intrinsic camera parameter matrix
%
% OUTPUTS:   pin       : [u;v] coordinates of the points used for H_in
%            pface     : [u;v] coordinates of the points used for H_f
%            pout      : [u;v] coordinates of the points used for H_out
%            ed1       : List containing the visible edges
%            see_face  : Vector containing the face-numbers which are
%                        visible
%% Determine points for H_in
%Rotate/translate the face centers and outer points to current state
outer_points = box.face_centers + (sign(box.face_centers)*0.01);
R_face_centers      = hom_trans(box.face_centers,state);
R_face_outer_points = hom_trans(outer_points,state);

%Face normal vectors and vectors from face center to the viewing point
normals     = R_face_outer_points-R_face_centers; 
view_vector = -R_face_centers; 

%Compute for each face the angle (theta) between normal and the vector to 
%the view-point. If the cos(theta) of that face is positive, the face is
%visible from the viewpoint.
visible=dot(normals,view_vector)./(vecnorm(normals).*vecnorm(view_vector));
see_face = find(visible>0);

%For each visible face, transform the points to the current state. Next, by
%perspective projection, compute their representation in the 2D plane.
pin = cell(1,6); pface = cell(1,6);
for ii=1:length(see_face)
pin{see_face(ii)}=persProj(hom_trans(box.Epoints{see_face(ii)},state),K);
pface{see_face(ii)}=persProj(hom_trans(box.Spoints{see_face(ii)},state),K);
end

%% Determine points for H_out
%Determine from the visible faces which edges are visible
model_faces{1}=[1,2,3,4];       model_faces{4}=[4,8,9,12];
model_faces{2}=[2,6,10,11];     model_faces{5}=[3,7,11,12];
model_faces{3}=[5,6,7,8];       model_faces{6}=[1,5,9,10];

%obtain the visible edges
ed = horzcat(model_faces{see_face});

%Only select the outside edges 
[counts, values] = histcounts(ed,1:13); %Count occurence of the edges
repeatedElements = values(counts >= 2); %Select edges that occure >= 2 
ed(ismember(ed,repeatedElements))=[];   %Remove these edges from the list

%Transform points of visible edges to the current state. Next, 
%by perspective projection, compute their representation in the 2D plane.
pts_out = hom_trans(horzcat(box.Opoints{ed}),state); 
pout = persProj(pts_out,K); 