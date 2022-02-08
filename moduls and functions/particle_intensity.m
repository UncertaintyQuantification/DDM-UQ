function [index_fill, Ip_fill]  = particle_intensity(xp,yp,sigma,Ic_here,L)

%xp=pos(1, 1,1);
%yp=pos(1, 1,2);

x_range=floor(xp-3*sigma):ceil(xp+3*sigma);
y_range=floor(yp-3*sigma):ceil(yp+3*sigma);

[x_I,y_J]=ndgrid(x_range,y_range);

binary_matrix=logical(((x_I-xp).^2+(y_J-yp).^2)<=((3*sigma)^2));

dist_2=((x_I-xp).^2+(y_J-yp).^2);

Ip=Ic_here * exp(-dist_2 / (2*sigma^2));

x_fill=x_I(binary_matrix);
y_fill=y_J(binary_matrix);

index_fill=x_fill+L*(y_fill-1);

Ip_fill=Ip(binary_matrix);
end
