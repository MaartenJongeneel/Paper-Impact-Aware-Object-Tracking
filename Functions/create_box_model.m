function [box] = create_box_model(l,w,h,m,doPlot)
%% Cuboid model
% This script computes the model for the cuboid.
%
% INPUTS:    l                 : Length of the box: X-direction [m]
%            w                 : Width of the box: Y-direction  [m]
%            h                 : Height of the box: Z-direction [m]
%            m                 : Mass ofthe box                 [kg]
%            doPlot            : Boolean deciding if you want to plot the
%                                geometric model of the box
%
% OUTPUTS:   box               : struct with fields:
%               dim            : Array containing dimensions of the box [l,w,h]
%               face_centers   : [x;y;z] coordinates of the face centers
%               Epoints        : Cell array with points for H_in for each face
%               Spoints        : Cell array with points for H_face for each face
%               Opoints        : Cell array with points for H_out for each edge
%               model_faces    : Array indicating relation between faces and edges
%% Creating the model
%%
% Constants and settings
pde    = 0.003;  %Distance point on surface from edge [m]
pd     = 0.005;  %Distance between points on surface  [m]
fd     = 0.015;  %Distance points from edge to surface points [m]
sd     = 0.003;  %Difference in size between the actual 3D model and the bigger one [m]
doSave = false;  %Decide if you want to save results  

%% Create the box struct
%Mass matrix of the box
Ml = m*eye(3);

%Inertia matrix of the box
I  = [(m/12)*(w^2+h^2),                 0,                  0;
                     0,  (m/12)*(l^2+h^2),                  0;
                     0,                 0,   (m/12)*(l^2+w^2);];
%Inertia tensor
box.B_M_B = [Ml zeros(3,3); zeros(3,3) I];

%Mass of the box
box.mass = m;

%% Define vertices for normal and bigger model and face centers
box.dim = [l,w,h];
box.vertices  = [l -l -l l l -l -l l; -w -w w w -w -w w w; -h -h -h -h h h h h]/2;
vert_b = box.vertices+sign(box.vertices)*sd;
box.face_centers = [0 -l 0 l 0 0; 0 0 0 0 w -w; -h 0 h 0 0 0]/2;

%% List the edges and faces of the model
%Define for each edge the corresponding vertices
edges = [1 2 3 4 5 6 7 8 1 2 3 4; 2 3 4 1 6 7 8 5 5 6 7 8];
%Define for each face the min/max values of the changing coordinates and the constant coordinate (5th column)
minmax=[-l  l  -w  w -h; -w  w  -h  h -l; -l  l  -w  w  h; -w  w  -h  h  l; -l  l  -h  h  w; -l  l  -h  h -w]/2;
    
%Index1 gives the corresponding direction for the values above(1=x,2=y,3=z)
index1 = [1,2,3;2,3,1;1,2,3;2,3,1;1,3,2;1,3,2];

%%  Define the 3D points where to sample color, for each face. 
for fi=1:length(box.face_centers) %For each face
    dis1=(minmax(fi,2)-pde)-(minmax(fi,1)+pde); %Compute distance between inner vertices
    dis2=(minmax(fi,4)-pde)-(minmax(fi,3)+pde); %Compute distance between inner vertices
    dis3=(minmax(fi,2)-fd)-(minmax(fi,1)+fd); %Compute distance between inner vertices
    dis4=(minmax(fi,4)-fd)-(minmax(fi,3)+fd); %Compute distance between inner vertices
    
    %Below, l1 and l2 are the linear spaces from which points are sampled on each face
    l1 = linspace((minmax(fi,1)+pde),(minmax(fi,2)-pde),round((dis1/pd))+1); 
    l2 = linspace((minmax(fi,3)+pde),(minmax(fi,4)-pde),round((dis2/pd))+1);
    l3 = linspace((minmax(fi,1)+fd),(minmax(fi,2)-fd),round((dis3/pd))+1); 
    l4 = linspace((minmax(fi,3)+fd),(minmax(fi,4)-fd),round((dis4/pd))+1);
    
    %Create sampled points near the edges
    [X,Y]=meshgrid(l1,l2);
    X(2:end-1,2:end-1) = NaN;
    Y(2:end-1,2:end-1) = NaN;
    pnts_e = [reshape(X,[],1),reshape(Y,[],1),repmat(minmax(fi,5),length(reshape(Y,[],1)),1)];
    pnts_e(any(isnan(pnts_e),2),:) = [];
    pnts_e(:,index1(fi,:)) = pnts_e;
    box.Epoints{fi}=pnts_e';
    
    %Create sampled points on the surface
    [X,Y]=meshgrid(l3,l4);
    X(2:end-1,2:end-1) = NaN;
    Y(2:end-1,2:end-1) = NaN;
    pnts_s = [reshape(X,[],1),reshape(Y,[],1),repmat(minmax(fi,5),length(reshape(Y,[],1)),1)];
    pnts_s(any(isnan(pnts_s),2),:) = [];
    pnts_s(:,index1(fi,:)) = pnts_s;
    box.Spoints{fi}=pnts_s';
end

%Compute the outer points 
for ei=1:length(edges)
    vec=vert_b(:,edges(1,ei))-vert_b(:,edges(2,ei));  %Difference between two vertices
    vecn=vec/norm(vec);                               %Normalized vector pointing from one vertex to another
    lin = linspace(0,norm(vec),((norm(vec)/pd)+1));   %Create lin space: points between the vertices
    box.Opoints{ei}=vert_b(:,edges(2,ei))+vecn*lin;    %Store them in a cell aray
end

%% Save the data
if doSave
    save('box.mat','box');
end

%% Plotting Figure
if doPlot
    set(groot,'defaulttextinterpreter','latex'); set(groot,'defaultAxesTickLabelInterpreter','latex');  set(groot,'defaultLegendInterpreter','latex');
    figure;
    hold on;
    for fi = 1:6
        plot3(box.Epoints{fi}(1,:), box.Epoints{fi}(2,:), box.Epoints{fi}(3,:),'.','color',[0 0.4470 0.7410]);
        plot3(box.Spoints{fi}(1,:), box.Spoints{fi}(2,:), box.Spoints{fi}(3,:),'.','color',[0.4660 0.6740 0.1880]);
    end
    for jj=1:12
        plot3(box.Opoints{jj}(1,:),box.Opoints{jj}(2,:),box.Opoints{jj}(3,:),'+','color',[0.8500 0.3250 0.0980],'linewidth',2);
        %Plot the black contour of the box
        two_vertices=[box.vertices(:,edges(1,jj)),box.vertices(:,edges(2,jj))];
        line(two_vertices(1,:),two_vertices(2,:),two_vertices(3,:),'color','k','LineWidth',2);
    end
    axis equal
    view(-38,34)
    axis([-0.060 0.060 -0.085 0.085 -0.035 0.035]);
    xlabel('x [m]','position',[  0 -0.100 -0.035]);
    ylabel('y [m]','position',[-0.075    0 -0.035]);
    zlabel('z [m]','position',[-0.070  0.100   0]);
    grid on;
end
  