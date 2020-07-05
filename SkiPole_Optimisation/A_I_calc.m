%% Engineering optimization 
% Area and intertia calculation
% Martin Janssens
% Rens Liebrand

function [A,I_x,I_y,x_bar,y_bar]=A_I_calc(r,h,c,t)
% Import coordinates and close the airfoil, just blund
[x_coord,y_coord]=Airfoil_coordinates(r,h,c);
x_coord=[x_coord; x_coord(1)];
y_coord=[y_coord; y_coord(1)];

% Determine element lengths
dx=abs(x_coord(2:end)-x_coord(1:end-1));
dy=abs(y_coord(2:end)-y_coord(1:end-1));
element_length=sqrt(dx.^2+dy.^2);

% Assuming constant t and no overlap deterimne A
dA=element_length*t;
A=sum(dA);

% Determine centroid by xbar=(sum x*dA)/(sum dA)
x_element=(x_coord(2:end)+x_coord(1:end-1))/2;
y_element=(y_coord(2:end)+y_coord(1:end-1))/2;
x_bar=sum(x_element.*dA)/A;
y_bar=sum(y_element.*dA)/A;

% Determine moment of inertia 
% I0=Ibar_ + Ad^2
theta=90-atand(dy./dx); % See drawing, is defined the other way around
Ibar_x=t*element_length./12.*(element_length.^2.*cosd(theta)+t^2.*sind(theta));
Ibar_y=t*element_length./12.*(element_length.^2.*sind(theta)+t^2.*cosd(theta));
% Add Steiner terms
dx_centroid=abs(abs(x_element)-x_bar);
dy_centroid=abs(abs(y_element)-y_bar);
I_x=sum(Ibar_x+dA.*dx_centroid.^2);
I_y=sum(Ibar_y+dA.*dx_centroid.^2);

%% Plot Panels
% figure()
% hold on
% plot(x_coord*1000,y_coord*1000,'k-','Linewidth',5)
% scatter(x_coord*1000,y_coord*1000,'g')
% scatter(x_bar*1000,y_bar*1000)
% xlabel('$x$ [mm]','fontsize',16,'Interpreter','LaTex')
% ylabel('$y$ [mm]','fontsize',16,'Interpreter','LaTex')
% leg=legend('Lineairised Airfoil Shape','Coordinates','Centroid','location','northeast');
% leg.FontSize = 12;
end
