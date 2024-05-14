function [] = plot_contourtop_drift(Contopdata,plot_type,FIG,savename)

% if exist(savename.pcb,"file")
%     load(savename.pcb,'PCBgrid2D','PCBdata2D')
% else
%     warning([savename.pcb, 'does not exist.'])
% end
if exist(savename.highmesh_psi,"file")
    load(savename.highmesh_psi,'highmesh_PCBgrid2D','highmesh_PCBdata2D')
else
    warning([savename.highmesh_psi, 'does not exist.'])
end

figure('Position',[0 0 1500 1500],'visible','on')
for i_t = 1:FIG.tate*FIG.yoko
    subplot(FIG.tate,round(FIG.yoko),i_t)
    switch plot_type
        case 'V_r'
            plot(Contopdata.rq,Contopdata.Vcontop_r(:,i_t))
            xlabel('R [m]')
            ylabel('Psi Contour Flow [km/s]')
        case 'scatter'
            if exist(savename.pcb,"file")
                hold on
                idx_pcb = knnsearch(highmesh_PCBdata2D.trange',FIG.start+(i_t-1)*FIG.dt);
                contour(highmesh_PCBgrid2D.zq(1,:),highmesh_PCBgrid2D.rq(:,1),squeeze(highmesh_PCBdata2D.psi(:,:,idx_pcb)),[-20e-3:0.1e-3:40e-3],'black','LineWidth',1)
            end
            scatter(Contopdata.zq(:,i_t),Contopdata.rq,[],Contopdata.Vcontop_r(:,i_t),"filled")
            daspect([1 1 1])
            xlim([-0.02 0.05])
            % ylim([0.1 0.27])
            xlabel('Z [m]')
            ylabel('R [m]')
            view([90 -90])%RZ反転
            colorbar
            clim([-30 20])
    end
    title([num2str(Contopdata.trange(i_t)) 'us'])
    ax = gca;
    ax.FontSize = 12;
end

end