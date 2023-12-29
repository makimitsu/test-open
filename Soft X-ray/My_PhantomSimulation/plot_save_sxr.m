function plot_save_sxr(range,EE,EEorigin)

range = range./1000;
zmin1 = range(1);
zmax1 = range(2);
zmin2 = range(3);
zmax2 = range(4);
rmin = range(5);
rmax = range(6);
r_space_SXR = linspace(rmin,rmax,size(EE,1));
z_space_SXR1 = linspace(zmin1,zmax1,size(EE,2));
z_space_SXR2 = linspace(zmin2,zmax2,size(EE,2));

r_range = find(0.060<=r_space_SXR & r_space_SXR<=0.330);
r_space_SXR = r_space_SXR(r_range);
z_range1 = find(-0.12<=z_space_SXR1 & z_space_SXR1<=0.12);
z_range2 = find(-0.12<=z_space_SXR2 & z_space_SXR2<=0.12);

z_space_SXR1 = z_space_SXR1(z_range1);
z_space_SXR2 = z_space_SXR2(z_range2);

% if show_flux_surface
%     psi_mesh_z = grid2D.zq;
%     psi_mesh_r = grid2D.rq;
%     t_idx = find(data2D.trange==t);
%     psi = data2D.psi(:,:,t_idx);
%     
%     psi_min = min(min(psi));
%     psi_max = max(max(psi));
%     contour_layer = linspace(psi_min,psi_max,20);
% end

[SXR_mesh_z1,SXR_mesh_r] = meshgrid(z_space_SXR1,r_space_SXR);
[SXR_mesh_z2,~] = meshgrid(z_space_SXR2,r_space_SXR);

% 負の要素を0で置換
negativeEE = find(EE<0);
EE(negativeEE) = zeros(size(negativeEE));

% if i <= 2
    % z_range = z_range2;
    % SXR_mesh_z = SXR_mesh_z2;
% else
    z_range = z_range1;
    SXR_mesh_z = SXR_mesh_z1;
% end
EE_plot = EE(r_range,z_range);
figure(1);
[~,h] = contourf(SXR_mesh_z,SXR_mesh_r,EE_plot,20);colormap("parula");axis([-0.12 0.12 0.06 0.33]);daspect([1 0.8 1]);
h.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=12;clim([0 0.04]);
set(gcf,'Name','再構成結果','NumberTitle','off');

EEorigin_plot = EEorigin(r_range,z_range);
figure(2);
[~,h] = contourf(SXR_mesh_z,SXR_mesh_r,EEorigin_plot,20);colormap("parula");axis([-0.12 0.12 0.06 0.33]);daspect([1 0.8 1]);
h.LineStyle = 'none';
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=12;clim([0 0.05]);
set(gcf,'Name','元画像','NumberTitle','off');
end