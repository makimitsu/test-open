function [] = movie_ESP(plot_Efield,trange,data2D)
figure('Position', [0 0 1500 1500],'visible','on');
[~,h] = contourf(data2D.phi_mesh_z,data2D.phi_mesh_r,squeeze(data2D.phi_grid(1,:,:)),'edgecolor','none');
colormap(redblue(300));
colorbar
xlim([-0.2 0.2])
ylim([0.1 0.3])
clim([-200 200])
hold on
if plot_Efield
    q = quiver(data2D.phi_mesh_z,data2D.phi_mesh_r,squeeze(data2D.Ez_grid(1,:,:)),squeeze(data2D.Er_grid(1,:,:)),3);
    q.Color = "k";
end
title([num2str(trange(1)) 'us'])
for i = 2:numel(trange)
    if mod(i-1,5) == 0
        h.ZData = squeeze(data2D.phi_grid(i,:,:));
        if plot_Efield
            q.UData = squeeze(data2D.Ez_grid(i,:,:));
            q.VData = squeeze(data2D.Er_grid(i,:,:));
        end
        title([num2str(trange(i)) 'us'])
        drawnow
    end
end
end

