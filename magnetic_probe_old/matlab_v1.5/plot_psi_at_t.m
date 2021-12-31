function [] = plot_psi_at_t(B_z,r_probe,z_probe,t,fitting,fill,fixed_Clayer,show_probe)
% plot psi in rz plane
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
%   integer: t, time of interest (us)
%   boolean: fitting, option for enabling cubic spline fitting
%   boolean: fixed_Clayer, option for enabling 
%            true -> fixed value between each contour layer
%            false-> fixed number of total layers = layer_num
%   boolean: show_probe, option for showing location of pickup coils

layer_num = 40;
% a uniform mesh for interpolation
[probe_mesh_z,probe_mesh_r] = meshgrid(z_probe,r_probe);

if fitting
    z_space = linspace(z_probe(1),z_probe(end),50);
    r_space = linspace(r_probe(1),r_probe(end),50);
    [psi_mesh_z,psi_mesh_r] = meshgrid(z_space,r_space);
else
    [psi_mesh_z,psi_mesh_r] = meshgrid(z_probe,r_probe);
end

if fixed_Clayer
    minx = round(min(B_z(:)),3);
    maxx = round(max(B_z(:)),3);
    contour_layer =  minx:0.001:maxx;
else
    contour_layer = layer_num;
end

figure('Position', [0 0 1200 600])

hold on
psi = get_psi(B_z,r_probe,t);

if fitting
    psi = griddata(z_probe,r_probe,psi,psi_mesh_z,psi_mesh_r,'cubic');
end

if fill
    contourf(psi_mesh_z,psi_mesh_r,psi,contour_layer,'Fill','on');
else
    contourf(psi_mesh_z,psi_mesh_r,psi,contour_layer,'Fill','off');
end

if show_probe
    plot(probe_mesh_z,probe_mesh_r,'*','color','k','markersize',3)
end
    
title(strcat(num2str(t),' us','; Psi (Wb)'));
xlabel('z (m)');
ylabel('r (m)');

hold off
end