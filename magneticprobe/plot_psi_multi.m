function [] = plot_psi_multi(B_z,r_probe,z_probe,times,fitting,fill,fixed_Clayer,show_probe,shot)
% plot psi in rz plane
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
%   1d array of integer: t, times of interest (us) e.g.)[461:2:480]
%   boolean: fitting, option for enabling cubic spline fitting 
%   to fit the the 8x28 mesh to 50x50 mesh
%   boolean: fixed_Clayer, option for enabling 
%            true -> fixed value between each contour layer
%            false-> fixed number of total layers = layer_num
%   boolean: adjust, option for showing location of pickup coils

layer_num = 20;
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
    max_psi_color = 0.0400;
    min_psi_color = -0.0200;
    layer_resolution = 0.0002;
    contour_layer =  min_psi_color:layer_resolution:max_psi_color;
else
    contour_layer = layer_num;
end

maxRow = 5;
cnt = 1; 
column = 5; 
row = ceil(numel(times)/column);
if row > maxRow
    row = maxRow;
end
f = figure('name',['shot', num2str(shot)]);
f.WindowState = 'maximized';
for t = times
    if cnt > maxRow*column
        f = figure;
        f.WindowState = 'maximized';
        cnt = 1;
    end
    subplot(row,column,cnt)
    %figure('Position', [0 0 1200 600])
    hold on

    psi = get_psi(B_z,r_probe,t);

    if fitting
        psi = griddata(z_probe,r_probe,psi,psi_mesh_z,psi_mesh_r,'cubic');
    end

    if fill
        contourf(psi_mesh_z,psi_mesh_r,psi,contour_layer,'Fill','on');
    else
        contourf(psi_mesh_z,psi_mesh_r,psi,contour_layer,'-k','Fill','off');
    end

    if show_probe
        plot(probe_mesh_z,probe_mesh_r,'*','color','k','markersize',3)
    end
    

    title(strcat(num2str(t),' us','; Psi (Wb)'));
    xlabel('z (m)');
    ylabel('r (m)');
    colorbar;
    ax = gca;
    ax.FontSize = 5; 
    hold off
    daspect([1 1 1])
    cnt = cnt + 1;
end
end