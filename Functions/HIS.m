function [His] = HIS(img,Bh,Bs,Bi)
% This function computes the color histogram of an reference image
%
% INPUTS:    img  :   Image data obtained by imread('file')
%            Bh   :   Number of bins used for hue histogram
%            Bs   :   Number of bins used for saturation histogram
%            Bi   :   Number of bins used for intensity histogram
%
% OUTPUTS:   His  :   Normalized color histogram with dimensions Bh*Bs*Bi

%% Obtaining the histogram
[height,width] = size(img(:,:,1));   %Height and width of the image
hsi = rgb2hsi(img,Bh,Bs,Bi);         %Obtain HSI values of image

His = zeros(12,12,4);                %Preallocate histogram
for x = 1:width
    for y = 1:height
        if (round((double(img(y,x,1))+double(img(y,x,2))+double(img(y,x,3)))/3)<249) %White pixels are not used
            H=hsi(y,x,1);
            S=hsi(y,x,2);
            I=hsi(y,x,3);
            His(H,S,I)=His(H,S,I)+1;
        end
    end
end

total=sum(sum(sum(His)));
His=His./total;   %Normalize the histogram (range 0 to 1)
