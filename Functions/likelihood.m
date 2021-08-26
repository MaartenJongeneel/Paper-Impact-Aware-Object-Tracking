function [L,S1,S2,S3]=likelihood(state,Href,hsi,K,box)
% Computes the likelihood 
% INPUTS:    state        : cell(R;o;v;w): state of the particle 
%            Href         : Cell aray with reference histogram of each face
%            hsi          : HSI image of the current frame
%            K            : Intrinsic camera matrix (3x3)
%            box          : Struct containing geometric properties of the
%                           box
%
% OUTPUTS:   L            : Likelihood 
%            S1           : Similarity 1
%            S2           : Similarity 2
%            S3           : Similarity 3
%% Compute the likelihood
%Compute the points on the surface of the cuboid and outside the cuboid
[pin,pface,pout,see_face] = compute_points(box,state,K);
N = length(see_face);

%Compute histograms for both sets of points
for ii = 1:N
Hin{see_face(ii)}  = Hpoints(pin{see_face(ii)},hsi); 
Hf{see_face(ii)} = Hpoints(pface{see_face(ii)},hsi);
end

%Weighting parameters for the likelihood function
k1=1; k2=2; k3=5; b=0.1;

[Hout] = Hpoints(pout,hsi); 

%Compute the likelihood based on Bhattacharyya distance between histograms.
S1(1:3)=0;
S2(1:3)=0;
S3(1:3)=0;
for ii = 1:N
    for c1=1:12
        for c2=1:12
            for c3=1:4
                S1(ii)=S1(ii)+ sqrt(Hin{see_face(ii)}(c1,c2,c3)*Href{see_face(ii)}(c1,c2,c3)); %Similarity contour inside points and reference
                S2(ii)=S2(ii)+ sqrt(Hf{see_face(ii)}(c1,c2,c3)*Href{see_face(ii)}(c1,c2,c3));  %Similarity face points and reference
                S3(ii)=S3(ii)+ sqrt(Hin{see_face(ii)}(c1,c2,c3)*Hout(c1,c2,c3));               %Similarity contour inside and contour outside points               
            end
        end
    end
end
S1 = sum(S1(:))/N;
S2 = sum(S2(:))/N;
S3 = sum(S3(:))/N;
%Specify the distance function
D = (k1*(1-S1)+k2*(1-S2)+k3*S3)/(k1+k2+k3); 
%Likelihood function
L = exp(-(abs(D))/b);
