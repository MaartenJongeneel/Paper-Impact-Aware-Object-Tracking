function [His]=Hpoints(points,hsi)
% Compute the color histogram of the points in the image hsi
% INPUTS:    points:  [x;y] array of points from which the histogram needs
%                     to be computed
%            hsi   :  HSI image with Hue Saturation and Intensity bins of
%                     12 12 and 4 respectively
% OUTPUTS:   His   :  12x12x4 Color Histogram
%% Script
[~,Npoints]      = size(points);              %Compute the number of points
[rows,columns,~] = size(hsi);                 %Obtain the size of the image
His=zeros(12,12,4);                           %Preallocate memory for speed

count=0;
for ii=1:Npoints
    u=round(points(1,ii));                    %x position in the image
    v=round(points(2,ii));                    %y position in the image
    if((v<=rows)&&(u<=columns)&&(v>=1)&&(u>=1))
        H=hsi(v,u,1);                         %H value at position v,u
        S=hsi(v,u,2);                         %S value at position v,u
        I=hsi(v,u,3);                         %I value at position v,u
        His(H,S,I)=His(H,S,I)+1;              %Store the value
        count=count+1;
    end
end

%Normalize the histogram by dividing each value by the total sum
total=sum(sum(sum(His)));
if(total~=0)
    His=His./total;                 
end
