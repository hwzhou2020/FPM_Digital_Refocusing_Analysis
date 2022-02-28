function [ uSobel ] = Sobel( image )
%SOBEL �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[high,width] = size(image);   % ���ͼ��ĸ߶ȺͿ��   
U = double(image);   
uSobel = image;
for i = 2:high - 1   %sobel��Ե���
    for j = 2:width - 1
        Gx = (U(i+1,j-1) + 2*U(i+1,j) + U(i+1,j+1)) - (U(i-1,j-1) + 2*U(i-1,j) + U(i-1,j+1));
        Gy = (U(i-1,j+1) + 2*U(i,j+1) + U(i+1,j+1)) - (U(i-1,j-1) + 2*U(i,j-1) + U(i+1,j-1));
        uSobel(i,j) = sqrt(Gx^2 + Gy^2); 
    end
end 
end

