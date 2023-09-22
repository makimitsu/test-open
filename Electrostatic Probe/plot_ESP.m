function [] = plot_ESP(ESP,data2D)
figure('Position', [0 0 1500 1500],'visible','on');
for i = 1:ESP.tate*ESP.yoko
    offset_t = knnsearch(ESP.trange',ESP.start_t);
    idx_t = offset_t+(i-1)*ESP.dt*10;
    subplot(ESP.tate,ESP.yoko,i)
    contourf(data2D.phi_mesh_z,data2D.phi_mesh_r,squeeze(data2D.phi_grid(idx_t,:,:)),'edgecolor','none');
    colormap(redblue(3000));
    c = colorbar;
    % xlim([-0.2 0.2])
    % ylim([0.1 0.3])
    xlim([-0.1 0.1])
    % xlim([-0.17 0.17])
    ylim([0.07 0.3])
    pbaspect([1 1 1])
    clim([-240 240])
    c.Label.String = 'Floating phi [V]';
    hold on
    if ESP.vector
        q = quiver(data2D.phi_mesh_z,data2D.phi_mesh_r,squeeze(data2D.Ez_grid(idx_t,:,:)),squeeze(data2D.Er_grid(idx_t,:,:)),3);
        q.Color = "k";
    end
    hold on
    title([num2str(ESP.trange(idx_t)) 'us'])
end
end

