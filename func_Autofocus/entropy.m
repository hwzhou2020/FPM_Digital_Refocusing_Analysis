function [k] = entropy( X )
%ENTROPY �����άͼ�����
% X ����Ķ�άͼ��
% k �������ֵ

X = im2uint8(mat2gray(X));
[C,R]=size(X);      %��ͼ��Ĺ��
Img_size=C*R;    %ͼ�����ص���ܸ���
L=256;                %ͼ��ĻҶȼ�
H_img=0;
nk=zeros(L,1);
for i=1:C
    for j=1:R
        Img_level=X(i,j)+1;                 %��ȡͼ��ĻҶȼ�           
        nk(Img_level)=nk(Img_level)+1;      %ͳ��ÿ���Ҷȼ����صĵ���
    end
end
for k=1:L
    Ps(k)=nk(k)/Img_size;                  %����ÿһ���Ҷȼ����ص���ռ�ĸ���
    if Ps(k)~=0                           %ȥ������Ϊ0�����ص�
    H_img=-Ps(k)*log2(Ps(k))+H_img;        %����ֵ�Ĺ�ʽ
    end
end
k=H_img;
end

