function [] = plot_magline_drift(Maglinedata2D,plot_type,FIG,savename)

% if exist(savename.highmesh_psi,"file")
%     load(savename.highmesh_psi,'highmesh_PCBgrid2D','highmesh_PCBdata2D')
% else
%     warning([savename.highmesh_psi, 'does not exist.'])
% end
if exist(savename.pcb,"file")
    load(savename.pcb,'PCBgrid2D','PCBdata2D')
else
    warning([savename.pcb, 'does not exist.'])
end

figure('Position',[0 0 1500 1500],'visible','on')
for i_t = 1:FIG.tate*FIG.yoko
    subplot(FIG.tate,round(FIG.yoko),i_t)
    switch plot_type
        case 'V_r'
            contourf(Maglinedata2D.zq,Maglinedata2D.rq,squeeze(Maglinedata2D.Vmagline_r(:,:,i_t)),100,'edgecolor','none')
            c = colorbar;
            colormap(jet)
            clim([min(Maglinedata2D.Vmagline_r,[],"all") max(Maglinedata2D.Vmagline_r,[],"all")])
            c.Label.String = 'Magnetic Line Drift [km/s]';
        case 'V_zr'
            q = quiver(Maglinedata2D.zq,Maglinedata2D.rq,squeeze(Maglinedata2D.Vmagline_z(:,:,i_t)),squeeze(Maglinedata2D.Vmagline_r(:,:,i_t)));
            q.Color = 'b';
            q.LineWidth = 2;
            q.AutoScaleFactor = 0.5;
    end
    % if exist(savename.highmesh_psi,"file")
    %     hold on
    %     idx_pcb = knnsearch(highmesh_PCBdata2D.trange',FIG.start+(i_t-1)*FIG.dt);
    %     contour(highmesh_PCBgrid2D.zq(1,:),highmesh_PCBgrid2D.rq(:,1),squeeze(highmesh_PCBdata2D.psi(:,:,idx_pcb)),[-20e-3:0.05e-3:40e-3],'black','LineWidth',1)
    % end
    if exist(savename.pcb,"file")
        hold on
        idx_pcb = knnsearch(PCBdata2D.trange',FIG.start+(i_t-1)*FIG.dt);
        contour(PCBgrid2D.zq(1,:),PCBgrid2D.rq(:,1),squeeze(PCBdata2D.psi(:,:,idx_pcb)),[-20e-3:0.1e-3:40e-3],'black','LineWidth',1)
    end
    title([num2str(Maglinedata2D.trange(i_t)) 'us'])
    daspect([1 1 1])
    % xlim([z_min z_max])
    % ylim([r_min r_max])
    xlim([-0.05 0.08])
    ylim([0.1 0.27])
    xlabel('Z [m]')
    ylabel('R [m]')
    ax = gca;
    % ax.FontSize = 60;
end

% switch plot_type
%     case {'F_r','F_z','F','V_r','V_z','V'}
%         sgtitle(c.Label.String)
%     case 'F_zr'
%         sgtitle('Centrifugal Force Vector')
%     case 'V_zr'
%         sgtitle('Curvature Drift Vector')
% end

% view([90 -90])%RZ反転
% clim([0 1.2E-17])
% clim([-1E-18 1E-18])
% xlim([-0.01 0.05])
% c.Location = "north";