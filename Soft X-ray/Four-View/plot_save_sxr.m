function plot_save_sxr(range,EE1,EE2,EE3,EE4)
% この関数は，発光出力EEをrangeで定義される範囲で切り取り，figureで出力します

range = range./1000; % 単位をmmからmに変換します
zmin = range(1);
zmax = range(2);
rmin = range(5);
rmax = range(6);
r_space_SXR = linspace(rmin,rmax,size(EE1,1));
z_space_SXR = linspace(zmin,zmax,size(EE1,2));

r_range = find(0.060<=r_space_SXR & r_space_SXR<=0.330);
r_space_SXR = r_space_SXR(r_range);
z_range = find(-0.17<=z_space_SXR & z_space_SXR<=0.17);

z_space_SXR = z_space_SXR(z_range);

EE1 = EE1(r_range,z_range);
EE2 = EE2(r_range,z_range);
EE3 = EE3(r_range,z_range);
EE4 = EE4(r_range,z_range);

subplot(2,2,1);
[SXR_mesh_z,SXR_mesh_r] = meshgrid(z_space_SXR,r_space_SXR);
[~,h1] = contourf(SXR_mesh_z,SXR_mesh_r,EE1,20);
h1.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('1');

subplot(2,2,2);
[~,h2] = contourf(SXR_mesh_z,SXR_mesh_r,EE2,20);
h2.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('2');

subplot(2,2,3);
[~,h3] = contourf(SXR_mesh_z,SXR_mesh_r,EE3,20);
h3.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('3');

subplot(2,2,4);
[~,h4] = contourf(SXR_mesh_z,SXR_mesh_r,EE4,20);
h4.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('4');

end