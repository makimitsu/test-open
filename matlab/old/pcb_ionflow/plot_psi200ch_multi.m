%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%プロット枚数、プロット開始時間、
%プロット時間間隔を指定して磁気面をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_psi200ch_multi(n_plot,t_start,dt,trange,cmap,cbar,grid2D,data2D,ok_z,ok_r)

fig = figure('Position', [0 0 1500 1500],'visible','on');
for m=1:n_plot %図示する時間
    i=(t_start-trange(1)+1)+(m-1)*dt; %end
    t=trange(i);
    if n_plot == 1
    elseif n_plot == 4
        subplot(2,2,m);
    elseif mod(n_plot,4) == 0
        subplot(n_plot/4,4,m)
    elseif n_plot < 4
        subplot(1,n_plot,m)
    else
        close(fig)
        warning('Input error in n_plot.')%n_plotの入力エラー
        return;
    end
    if cmap
        switch cmap
            case 'Bz'
                contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i)*1e3,30,'LineStyle','none')%Bz
                clim([-50,50])%Bz
                if cbar
                    colorbar('Location','eastoutside')
                    c = colorbar;
                    c.Label.String = 'Bz [mT]';
                end
            case 'Br'
                contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Br(:,:,i)*1e3,30,'LineStyle','none')%Br
                clim([-50,50])%Br
                if cbar
                    colorbar('Location','eastoutside')
                    c = colorbar;
                    c.Label.String = 'Br [mT]';
                end
            case 'psi'
                contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i)*1e3,40,'LineStyle','none')%psi
                clim([-5,5])%psi
                if cbar
                    colorbar('Location','eastoutside')
                    c = colorbar;
                    c.Label.String = 'Psi [Wb]';
                end
            case 'Bt'
                contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),-100e-3:0.5e-3:100e-3,'LineStyle','none')%Bt
                clim([-0.1,0.1])%Bt
                if cbar
                    colorbar('Location','eastoutside')
                    c = colorbar;
                    c.Label.String = 'Bt [T]';
                end
            case 'Jt'
                contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none')%jt
                clim([-0.8*1e+6,0.8*1e+6]) %jt
                if cbar
                    colorbar('Location','eastoutside')
                    c = colorbar;
                    c.Label.String = 'Jt [A/m^{2}]';
                end
            case 'Et'
                contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none')%Et
                clim([-500,400])%Et
                if cbar
                    colorbar('Location','eastoutside')
                    c = colorbar;
                    c.Label.String = 'Et [V/m]';
                end
            otherwise
                close(fig)
                warning('Input error in cmap.')%cmapの入力エラー
                return;
        end
        colormap(jet)
        axis image
        axis tight manual
        %     colorbar('Location','eastoutside')
        %カラーバーのラベル付け
        %     c = colorbar;
        %     c.Label.String = 'Jt [A/m^{2}]';
        hold on
    end
    %     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
    %     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
    %     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
    contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),-20e-3:0.2e-3:40e-3,'black','LineWidth',1)
    %     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
    %     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
    hold on
    plot(ok_z,ok_r,"k.",'MarkerSize', 6)%測定位置
    hold off
    title(string(t)+' us')
    xlabel('z [m]')
    ylabel('r [m]')
end
% figure('Position', [0 0 1500 1500],'visible','on');
% for m=1:16 %図示する時間
%     i=(t_start-trange(1)+1)+(m-1)*dt; %end
%     t=trange(i);
%     subplot(4,4,m)
%     plot(grid2D.zq(1,:),data2D.Bt(25,:,i))
%     title(string(t)+' us')
%     ylim([-0.03 0.03])
% end
% saveas(gcf,strcat('/Users/yunhancai/Downloads/files/a039_',num2str(shot)))
% saveas(gcf,strcat('/Users/yunhancai/Downloads/a039_',num2str(shot),'.png'))
% close

