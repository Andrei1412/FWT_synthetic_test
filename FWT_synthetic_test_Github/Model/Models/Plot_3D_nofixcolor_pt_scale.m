%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vsc_min=100;Vsc_max=800;
% nu=0.33;
% Vpc_min=Vsc_min*sqrt(2*(1-nu)/(1-2*nu));
% Vpc_max=Vsc_max*sqrt(2*(1-nu)/(1-2*nu));
 
fig=figure(80);
% imagesc(squeeze(Vs(:,:,8)));
xslice=(round(nx/2)-1)*dx;
zslice=(round(nz/2)-1)*dx;
yslice=(round(ny/2)-1)*dx;
c1=round(ny/2);

Vs_mean=mean(Vs(:));
Vp_mean=mean(Vp(:));

subplot(121)
Vs_r=rot_xzy2xyz(Vs);
slice(y,x,z,Vs_r,yslice,xslice,[]);
view(45,45);
grid on;
daspect([1 1 1]);
colormap(jet(64));
colorbar('vertical');
shading interp;
title('Vs [km/s]')
% title('Gra__Vs')
% set(gca,'YDir','reverse')
axis([min(y) max(y)  min(x) max(x) min(z) max(z)]);
ylabel('x-axis [km]')
xlabel('y-axis [km]')
zlabel('z-axis [km]')
set(gca,'ZDir','reverse')
set(gca,'CLim',[0.98*Vs_mean 1.02*Vs_mean])
 
subplot(122)
Vp_r=rot_xzy2xyz(Vp);
slice(y,x,z,Vp_r,yslice,xslice,[]);
view(45,45);
grid on;
daspect([1 1 1]);
colormap(jet(64));
colorbar('vertical');
shading interp;
title('Vp [km/s]')
% title('Gra__Vp')
% set(gca,'YDir','reverse')
axis([min(y) max(y)  min(x) max(x) min(z) max(z)]);
ylabel('x-axis [km]')
xlabel('y-axis [km]')
zlabel('z-axis [km]')
set(gca,'ZDir','reverse')
set(gca,'CLim',[0.98*Vp_mean 1.02*Vp_mean])
 
