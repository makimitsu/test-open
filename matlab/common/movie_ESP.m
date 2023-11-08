function [] = movie_ESP(ESP,ESPdata2D,colorplot)
switch colorplot
    case {'phi','Ez','Er'}
        figure('Position', [0 0 1500 1500],'visible','on');
        switch colorplot
            case 'phi'
                [~,h] = contourf(ESPdata2D.zq,ESPdata2D.rq,squeeze(ESPdata2D.phi(1,:,:)),100,'edgecolor','none');
                c = colorbar;
                clim([-240 240])
                c.Label.String = 'Floating phi [V]';
                colormap(redblue(3000))
                hold on
            case 'Ez'
                [~,h] = contourf(ESPdata2D.zq,ESPdata2D.rq,squeeze(ESPdata2D.Ez(1,:,:)),100,'edgecolor','none');
                c = colorbar;
                clim([-5000 5000])
                c.Label.String = 'E_z [V/m]';
                colormap(redblue(3000))
                hold on
            case 'Er'
                [~,h] = contourf(ESPdata2D.zq,ESPdata2D.rq,squeeze(ESPdata2D.Er(1,:,:)),100,'edgecolor','none');
                c = colorbar;
                clim([-5000 5000])
                c.Label.String = 'E_r [V/m]';
                colormap(redblue(3000))
                hold on
        end
        if ESP.vector
            q = quiver(ESPdata2D.zq,ESPdata2D.rq,squeeze(ESPdata2D.Ez(1,:,:)),squeeze(ESPdata2D.Er(1,:,:)),3);
            q.Color = "k";
        end
        hold on
        title([num2str(ESPdata2D.trange(1)) 'us'])
        xlim([-0.1 0.1])
        ylim([0.07 0.3])
        daspect([1 1 1])
    otherwise
        return
end

for i = 2:numel(ESP.trange)
    if mod(i-1,5) == 0
        switch colorplot
            case 'phi'
                h.ZData = squeeze(ESPdata2D.phi(i,:,:));
            case 'Ez'
                h.ZData = squeeze(ESPdata2D.Ez(i,:,:));
            case 'Er'
                h.ZData = squeeze(ESPdata2D.Er(i,:,:));
        end
        if ESP.vector
            q.UData = squeeze(ESPdata2D.Ez(i,:,:));
            q.VData = squeeze(ESPdata2D.Er(i,:,:));
        end
        title([num2str(ESP.trange(i)) 'us'])
        drawnow
    end
end
end

