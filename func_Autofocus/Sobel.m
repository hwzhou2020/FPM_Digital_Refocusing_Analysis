function [ uSobel ] = Sobel( image )
%SOBEL 此处显示有关此函数的摘要
%   此处显示详细说明
[high,width] = size(image);   % 获得图像的高度和宽度   
U = double(image);   
uSobel = image;
for i = 2:high - 1   %sobel边缘检测
    for j = 2:width - 1
        Gx = (U(i+1,j-1) + 2*U(i+1,j) + U(i+1,j+1)) - (U(i-1,j-1) + 2*U(i-1,j) + U(i-1,j+1));
        Gy = (U(i-1,j+1) + 2*U(i,j+1) + U(i+1,j+1)) - (U(i-1,j-1) + 2*U(i,j-1) + U(i+1,j-1));
        uSobel(i,j) = sqrt(Gx^2 + Gy^2); 
    end
end 
end

