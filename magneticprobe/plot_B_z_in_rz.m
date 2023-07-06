function [] = plot_B_z_in_rz(B_z,r_probe,z_probe)
% plot psi in rz plane
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
row = 4;
column = 4;
contour_layer = 30;
time = 420:4:480;
% a mesh of points corresponding to probe locations
[probe_mesh_z,probe_mesh_r] = meshgrid(z_probe,r_probe);

figure('Position', [0 0 1500 400])
n = 1;
for i = time
    subplot(row,column,n);
    
    contourf(probe_mesh_z,probe_mesh_r,B_z(:,:,i),contour_layer);
    
    title(strcat(num2str(i),' us','; B(T)'));
    xlabel('z (m)');
    ylabel('r (m)');
    colorbar;
    n = n+1;
end

end