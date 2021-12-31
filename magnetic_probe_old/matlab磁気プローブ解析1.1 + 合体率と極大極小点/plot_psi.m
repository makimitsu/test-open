function [] = plot_psi(B_z,r_probe,z_probe)
% plot psi in rz plane
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z

contour_layer = 30;
time = [470,472,474,476,478,480,482,484,486,488];

% a uniform mesh for interpolation
[psi_mesh_z,psi_mesh_r] = meshgrid(linspace(z_probe(25),z_probe(1),25),...
                                    linspace(r_probe(7),r_probe(1),25));

figure('Position', [0 0 1500 400])
n = 1;
for i = time
    psi = get_psi(B_z,r_probe,i);
    psi_interp = griddata(z_probe,r_probe,psi,psi_mesh_z,psi_mesh_r,'cubic');
    subplot(2,5,n);
    contourf(psi_mesh_z,psi_mesh_r,psi_interp,contour_layer);
    title(strcat(num2str(i),' us','; Psi (Wb)'));
    c = colorbar;
    xlabel('z (m)');
    ylabel('r (m)');
    colorbar;
    n = n+1;
end
end