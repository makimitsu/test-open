function [] = plot_ESP(tate,yoko,start_t,dt,plot_Efield,trange,data2D)
figure('Position', [0 0 1500 1500],'visible','on');
for i = 1:tate*yoko
    offset_t = knnsearch(trange',start_t);
    idx_t = offset_t+(i-1)*dt*10;
    subplot(tate,yoko,i)
    contourf(data2D.phi_mesh_z,data2D.phi_mesh_r,squeeze(data2D.phi_grid(idx_t,:,:)),'edgecolor','none');
    colormap(redblue(300));
    colorbar
    xlim([-0.2 0.2])
    ylim([0.1 0.3])
    clim([-200 200])
    hold on
    if plot_Efield
        q = quiver(data2D.phi_mesh_z,data2D.phi_mesh_r,squeeze(data2D.Ez_grid(idx_t,:,:)),squeeze(data2D.Er_grid(idx_t,:,:)),1);
        q.Color = "k";
    end
    hold off
    title([num2str(trange(idx_t)) 'us'])
end
end

