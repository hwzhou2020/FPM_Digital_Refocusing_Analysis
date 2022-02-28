
n0 = 1;
t1 = 2.02e-3;
n1 = 1.5997; % Polystyrene for lambda = 522.6e-9;
t2 = 500e-6;
n2 = 1.4115; % 1:10 PDMS for lambda = 522.6e-9???;
n3 = 1;

figure;
plot(sin_thetav,sin_thetah,'r*','MarkerSize',16);
axis equal;title('NA space');
hold on;

for k = 1:length(sin_thetav)

sin_theta0 = sin_thetav(k,1);
% n0*sin_theta0 = n1*sin_theta1 = n2*sin_theta2 = n3*sin_theta3
sin_theta1 = n0*sin_theta0/n1;
offset1 = t1*tan(asin(sin_theta1));
sin_theta2 = n0*sin_theta0/n2;
offset2 = t2*tan(asin(sin_theta2));
xoff = offset1+offset2;
sin_thetav_corr(k,1) = (hhled(k)*ds_led+img_center(1)+xoff)./dd(k);

sin_theta0 = sin_thetah(k,1);
sin_theta1 = n0*sin_theta0/n1;
offset1 = t1*tan(asin(sin_theta1));
sin_theta2 = n0*sin_theta0/n2;
offset2 = t2*tan(asin(sin_theta2));
yoff = offset1+offset2;
sin_thetah_corr(k,1) = (vvled(k)*ds_led+img_center(2)+yoff)./dd(k);

end
plot(sin_thetav_corr,sin_thetah_corr,'b*','MarkerSize',10);