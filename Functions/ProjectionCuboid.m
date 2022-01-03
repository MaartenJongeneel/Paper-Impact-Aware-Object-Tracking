function [figout] = ProjectionCuboid(img,K,box,state,y,z,s,width,height)
% This function creates an image from the cuboid in space. The background
% image is given by img, the model of the object by model such that the
% function plots the model at a given state in the image.
% 
% INPUTS:    img     : background image (grey image) 
%            K       : Intrinsic camera parameter matrix 
%            box     : model from create_cuboid_shape_model
%            state   : State of the box (rotation matrix and position vect)
%            y       : height of the contact surface             [m]
%            z       : distance to center of the contact surface [m]
%            s       : size of the contact surface (square)      [m]
%
%
% OUTPUTS:   figout  : image output (uint8)
%% Script
%Get the information out of the model
vertices       = box.vertices;

facecorners{1} = [vertices(:,1) vertices(:,2) vertices(:,3) vertices(:,4)];
facecorners{2} = [vertices(:,2) vertices(:,3) vertices(:,7) vertices(:,6)];
facecorners{3} = [vertices(:,5) vertices(:,6) vertices(:,7) vertices(:,8)];
facecorners{4} = [vertices(:,1) vertices(:,4) vertices(:,8) vertices(:,5)];
facecorners{5} = [vertices(:,3) vertices(:,4) vertices(:,8) vertices(:,7)];
facecorners{6} = [vertices(:,1) vertices(:,2) vertices(:,6) vertices(:,5)];

%% First part is to draw lines in the image
%Determine which edges are visible and which faces are visible
[~,~,~,see_face] = compute_points(box,state,K);

%% Second part is to plot the contact surface and the cuboid
%Get the 3D points of the contact surface
contactsurface = [[s;y;z-s],[-s;y;z-s],[-s;y;z+s],[s;y;z+s],[s;y;z-s]];
%Compute the 2D points of this surface
x1 = persProj(contactsurface,K);

%Define the position of the corners in the image plane
xyzpoints = horzcat(facecorners{see_face}); %3D points of the corners
xyzpoints = hom_trans(xyzpoints,state);     %Transpose to true position
uvpoints  = persProj(xyzpoints,K);
%Define the colors of each face
colorface = 255*[[1 0 0];[0 1 0];[0 1 1];[1 0 1];[1 1 0];[0 0 1]];

%Color the contact surface
q2 = repmat([1:width],1,height); %xvector of image coordinates
q3 = repelem([1:height],width);  %yvector of image coordinates
bs = inpolygon(q2,q3,x1(1,:),x1(2,:)); %Give pixels on the contact surface
Bs = find(bs==1);           %Find the coordinates of the contact surface
innx = q2(Bs);              %And get the x-coordinates...
inny = q3(Bs);              %and the y-coordinates.

%paint the contact surface white in the image
img=paint_pixelsUV(img,[innx;inny],[255,255,255]);

%Color the surfaces of the cuboid
if length(uvpoints)==4  %If only 1 face is visible
   face1 = inpolygon(q2,q3,uvpoints(1,1:4),uvpoints(2,1:4));
   face1 = find(face1==1);
   innx1 = q2(face1);
   inny1 = q3(face1);
   img=paint_pixelsUV(img,[innx1;inny1],colorface(see_face(1),:));
elseif length(uvpoints)==8 %If two faces are visible
   face1 = inpolygon(q2,q3,uvpoints(1,1:4),uvpoints(2,1:4));
   face2 = inpolygon(q2,q3,uvpoints(1,5:8),uvpoints(2,5:8));
   face1 = find(face1==1);
   face2 = find(face2==1);
   innx1 = q2(face1);
   inny1 = q3(face1);
   innx2 = q2(face2);
   inny2 = q3(face2);
   img=paint_pixelsUV(img,[innx1;inny1],colorface(see_face(1),:));
   img=paint_pixelsUV(img,[innx2;inny2],colorface(see_face(2),:));
elseif length(uvpoints) == 12 %If 3 faces (max) are visible
   face1 = inpolygon(q2,q3,uvpoints(1,1:4),uvpoints(2,1:4));
   face2 = inpolygon(q2,q3,uvpoints(1,5:8),uvpoints(2,5:8));
   face3 = inpolygon(q2,q3,uvpoints(1,9:12),uvpoints(2,9:12));
   face1 = find(face1==1);
   face2 = find(face2==1);
   face3 = find(face3==1);
   innx1 = q2(face1);
   inny1 = q3(face1);
   innx2 = q2(face2);
   inny2 = q3(face2);
   innx3 = q2(face3);
   inny3 = q3(face3);
   img=paint_pixelsUV(img,[innx1;inny1],colorface(see_face(1),:));
   img=paint_pixelsUV(img,[innx2;inny2],colorface(see_face(2),:));
   img=paint_pixelsUV(img,[innx3;inny3],colorface(see_face(3),:));
end 
figout = img;

function image=paint_pixelsUV(image,points,color)

[~,point_number]=size(points);
[rows, columns,~]=size(image);
for i=1:point_number
    u=round(points(1,i)); %horizontal.
    v=round(points(2,i)); %vertical.
    if((v<=rows)&&(u<=columns)&&(v>=1)&&(u>=1))
        image(v,u,:)=color;
    end
end
        