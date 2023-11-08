function plot_psi(PCBgrid2D,PCBdata2D,IDSP,FIG,colorplot)
figure('Position', [0 0 1500 1500],'visible','on')
for m=1:FIG.tate*FIG.yoko 
    i=FIG.start-PCBdata2D.trange(1)+1+(m-1)*FIG.dt;
    t=PCBdata2D.trange(i);
    subplot(FIG.tate,FIG.yoko,m)
    %カラープロット
    switch colorplot
        case 'psi'
            contourf(PCBgrid2D.zq(1,:),PCBgrid2D.rq(:,1),PCBdata2D.psi(:,:,i),80,'LineStyle','none')
            colormap(jet)
            clim([-10e-3,10e-3])%psi
            c = colorbar;
            c.Label.String = 'Psi [Wb]';
        case 'Bt'
            contourf(PCBgrid2D.zq(1,:),PCBgrid2D.rq(:,1),PCBdata2D.Bt(:,:,i),50,'LineStyle','none')
            colormap(jet)
            clim([0.1,0.4])%Bt
            c = colorbar;
            c.Label.String = 'B_t [T]';
        case 'Bt_ext'
            contourf(PCBgrid2D.zq(1,:),PCBgrid2D.rq(:,1),PCBdata2D.Bt_ext(:,:,i),50,'LineStyle','none')
            colormap(jet)
            clim([0.1,0.4])%Bt_ext
            c = colorbar;
            c.Label.String = 'B_t by TF cur. [T]';
        case 'Bt_plasma'
            contourf(PCBgrid2D.zq(1,:),PCBgrid2D.rq(:,1),PCBdata2D.Bt_plasma(:,:,i),50,'LineStyle','none')
            colormap(jet)
            clim([-0.01,0.04])%Bt_plasma
            c = colorbar;
            c.Label.String = 'B_t by plasma [T]';
        case 'Br'
            contourf(PCBgrid2D.zq(1,:),PCBgrid2D.rq(:,1),PCBdata2D.Br(:,:,i),100,'LineStyle','none')
            colormap(jet)
            clim([-0.04,0.04])%Br
            c = colorbar;
            c.Label.String = 'B_r [T]';
        case 'Bz'
            contourf(PCBgrid2D.zq(1,:),PCBgrid2D.rq(:,1),PCBdata2D.Bz(:,:,i),100,'LineStyle','none')
            colormap(jet)
            clim([-0.05,0.05])%Bz
            c = colorbar;
            c.Label.String = 'B_z [T]';
        case 'Et'
            contourf(PCBgrid2D.zq(1,:),PCBgrid2D.rq(:,1),PCBdata2D.Et(:,:,i),100,'LineStyle','none')
            colormap(jet)
            clim([-500,500])%Et
            c = colorbar;
            c.Label.String = 'E_t [V/m]';
        case 'Jt'
            contourf(PCBgrid2D.zq(1,:),PCBgrid2D.rq(:,1),PCBdata2D.Jt(:,:,i),30,'LineStyle','none')
            colormap(jet)
            clim([-1E6,1E6])%Jt
            c = colorbar;
            c.Label.String = 'Jt [A/m^{2}]';
    end
    hold on
    %磁気面
    contour(PCBgrid2D.zq(1,:),PCBgrid2D.rq(:,1),squeeze(PCBdata2D.psi(:,:,i)),[-20e-3:0.2e-3:40e-3],'black','LineWidth',1)
    % plot(PCBgrid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),PCBgrid2D.rq(opoint(:,:,i),1),"bo")
    % plot(PCBgrid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),PCBgrid2D.rq(xpoint(:,:,i),1),"bx")
    hold on
    for i_r = 1: size(PCBgrid2D.ok_bz_matrix,1)
        for i_z = 1: size(PCBgrid2D.ok_bz_matrix,2)
            % if PCBgrid2D.ok_bz_matrix(i_r,i_z) == 1
                p = plot(PCBgrid2D.zprobepcb(i_z),PCBgrid2D.rprobepcb(i_r),"b+");%測定位置
                p.LineWidth = 3;
                p.MarkerSize = 12;
            % end
        end
    end
    % [ok_z,ok_r] = meshgrid(PCBgrid2D.zprobepcb,PCBgrid2D.rprobepcb);
    % p = plot(ok_z,ok_r,"r+");%測定位置
    % p.LineWidth = 3;
    % p.MarkerSize = 12;
    % hold on
    % %IDSP計測点
    % plot(IDSP.z,IDSP.r,'r+')
    % hold on
    % title(string(t)+' us')
    xlim([-0.17 0.17])
    % ylim([0.07 0.3])
    daspect([1 1 1])
    xlabel('Z [m]')
    ylabel('R [m]')
    ax = gca;
    ax.FontSize = 20;
    % ax.XTickLabel = cell(size(ax.XTickLabel));
    % ax.YTickLabel = cell(size(ax.YTickLabel));
    view([90 -90])%RZ反転
end
end

