function [] = plot_ESP(ESP,ESPdata2D,colorplot)
switch colorplot
    case {'phi','Ez','Er'}
        figure('Position', [0 0 1500 1500],'visible','on');
        for i = 1:ESP.tate*ESP.yoko
            offset_t = knnsearch(ESPdata2D.trange',ESP.start);
            idx_t = offset_t+(i-1)*ESP.dt*10;
            subplot(ESP.tate,ESP.yoko,i)
            switch colorplot
                case 'phi'
                    contourf(ESPdata2D.zq,ESPdata2D.rq,squeeze(ESPdata2D.phi(idx_t,:,:)),100,'edgecolor','none');
                    c = colorbar;
                    clim([-240 240])
                    c.Label.String = 'Floating phi [V]';
                    colormap('bluewhitered');
                    hold on
                case 'Ez'
                    contourf(ESPdata2D.zq,ESPdata2D.rq,squeeze(ESPdata2D.Ez(idx_t,:,:)),100,'edgecolor','none');
                    c = colorbar;
                    clim([-5000 5000])
                    c.Label.String = 'E_z [V/m]';
                    colormap('bluewhitered');
                    hold on
                case 'Er'
                    contourf(ESPdata2D.zq,ESPdata2D.rq,squeeze(ESPdata2D.Er(idx_t,:,:)),100,'edgecolor','none');
                    c = colorbar;
                    clim([-5000 5000])
                    c.Label.String = 'E_r [V/m]';
                    colormap('bluewhitered');
                    hold on
            end
            if ESP.vector
                q = quiver(ESPdata2D.zq,ESPdata2D.rq,squeeze(ESPdata2D.Ez(idx_t,:,:)),squeeze(ESPdata2D.Er(idx_t,:,:)),3);
                q.Color = "k";
            end
            hold on
            title([num2str(ESPdata2D.trange(idx_t)) 'us'])
            xlim([-0.1 0.1])
            ylim([0.07 0.3])
            daspect([1 1 1])
        end
    otherwise
        return
end
end

