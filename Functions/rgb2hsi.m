function [h,s,i] = rgb2hsi(image,Bh,Bs,Bi)
%   RGB2HSI Convert red-green-blue colors to hue-saturation-intensity.
%   H = RGB2HSI(M) converts an RGB color map to an HSI color map.
%   Each map is a matrix with any number of rows, exactly three columns,
%   and elements in the interval specified by Bh,Bs,Bi.  The columns of the input matrix,
%   M, represent intensity of red, blue and green, respectively.  The
%   columns of the resulting output matrix, H, represent hue, saturation
%   and color value, respectively.
%
%   HSI = RGB2HSI(RGB) converts the RGB image RGB (3-D array) to the
%   equivalent HSI image HSI (3-D array).
%
%   CLASS SUPPORT
%   -------------
%   If the input is an RGB image, it can be of class uint8, uint16, or 
%   double; the output image is of class double.  If the input is a 
%   colormap, the input and output colormaps are both of class double.
% 
%   Undocumented syntaxes:
%
%   [H,S,I] = RGB2HSV(RGB) converts the RGB image RGB (3-D array) to
%   the equivalent HSI image where H,S and I are in range 0 to 255 
%
%   HSI = rgb2hsv(RGB,Bh,Bs,Bi) converts the RGB image into the equivalent
%   HSI image where H,S and I are are in range 0 to Bh, 0 to Bs and 0 to Bi
%   respectively
%
%   See Alvy Ray Smith, Color Gamut Transform Pairs, SIGGRAPH '78.
%   C. B. Moler, 8-17-86, 5-10-91, 2-2-92.
%      revised by C. Griffin for uint8 inputs 7-26-96
%      revised by P. Gravel for faster execution and less memory
%   Copyright 1984-2003 The MathWorks, Inc. 
%   $Revision: 5.15.4.1 $  $Date: 2003/05/01 20:43:53 $
%   
%   Adapted by Maarten Jongeneel on 2018/12/15


switch nargin
  case 1
     if isa(image, 'uint8')
        image = double(image) / 255; 
     elseif isa(image, 'uint16')
        image = double(image) / 65535;
     end
     Bh = 255;
     Bs = 255;
     Bi = 255;
  case 4
     if isa(image, 'uint8')
        image = double(image) / 255; 
     elseif isa(image, 'uint16')
        image = double(image) / 65535;
     end
    
   otherwise
      error('Wrong number of input arguments.');
end
  
threeD = (ndims(image)==3); % Determine if input includes a 3-D array
if threeD
  g = image(:,:,2); b = image(:,:,3); r = image(:,:,1);
  siz = size(r);
  r = r(:); g = g(:); b = b(:);
else
    error('Input image not a RGB image.')
end

i=(r+g+b)/3;
v = max(max(r,g),b);
s = zeros(size(v));
h = zeros(size(v));
s = (v - min(min(r,g),b));

z = ~s;
s = s + z;
k = find(r == v);
h(k) = (g(k) - b(k))./s(k);
k = find(g == v);
h(k) = 2 + (b(k) - r(k))./s(k);
k = find(b == v);
h(k) = 4 + (r(k) - g(k))./s(k);
h = h/6;
k = find(h < 0);
h(k) = h(k) + 1;
h=(~z).*h;

k = find(v);
s(k) = (~z(k)).*s(k)./v(k);
k = find(~v);
s(k) = 0;

if nargin <2
    h=h*255; %Values in range from 0 to 255
    s=s*255;
    i=i*255;
else
    h=round((Bh-1)*h)+1; %Values in range [1,2...,Bh-1,Bh]
    s=round((Bs-1)*s)+1; %values in range [1,2...,Bs-1,Bs]
    i=round((Bi-1)*i)+1; %values in range [1,2...,Bi-1,Bi]
end

if nargout<=1
  if (threeD || nargin==3)
    h = reshape(h,siz);
    s = reshape(s,siz);
    i = reshape(i,siz);
    h=cat(3,h,s,i);
  else
    h=[h s i];
  end
else
  h = reshape(h,siz);
  s = reshape(s,siz);
  i = reshape(i,siz);
end

