function MS8ND = AllDireDiff(Xf)
% mean of square of eight-neighbor difference 
% [I(x,y)-I(x+i,y+j)]^2,i=[-1,0,1],j=[-1,0,1]
% Moore neighborhood
[N,M]=size(Xf);
new_pic=zeros(N-2,M-2);
for i = 2:N-1
    for j = 2:M-1
    sum_row1=(Xf(i-1,j-1)-Xf(i,j))^2+(Xf(i-1,j)-Xf(i,j))^2+(Xf(i-1,j+1)-Xf(i,j))^2;
    sum_row2=(Xf(i,j-1)-Xf(i,j))^2+(Xf(i,j)-Xf(i,j))^2+(Xf(i,j+1)-Xf(i,j))^2;
    sum_row3=(Xf(i+1,j-1)-Xf(i,j))^2+(Xf(i+1,j)-Xf(i,j))^2+(Xf(i+1,j+1)-Xf(i,j))^2;
    new_pic(i-1,j-1)=sum_row1+sum_row2+sum_row3;
    end
end
MS8ND=mean2(new_pic);
end