function [] = plot_ionvdist(IDSPdata,IDSP,EXP,pathname,plot_type,save_fig)
%再構成速度分布関数を描画
[n_rows,~] = size(IDSPdata.ppoints);
time = IDSP.time;%プロット時間
if save_fig
    figure('Position',[0 0 1500 1500],'visible','off')
else
    figure('Position',[0 0 1500 1500],'visible','on')
end

sgtitle([num2str(IDSP.date),' shot', num2str(IDSP.shot),' ',num2str(time),'us (Z = ', num2str(IDSP.z(1,1)),' m)'])
for k = 1:IDSP.n_r
    subplot(2,4,k)
    draw_F = reshape(IDSPdata.F(:,k), [numel(IDSPdata.Vx),numel(IDSPdata.Vy)]);%描画用に整形
    draw_F = flipud(draw_F);
    [summit_y,summit_x] = find(draw_F == max(draw_F(:)));
    % fprintf('Summit is (Vz = %.1f, Vr = %.1f)[km/s].\n',IDSPdata.Vx(summit_Vx),IDSPdata.Vy(summit_Vy));
    levels = [draw_F(summit_y,summit_x)*0.1, draw_F(summit_y,summit_x)*0.2, ...
        draw_F(summit_y,summit_x)*0.3, draw_F(summit_y,summit_x)*0.4, ...
        draw_F(summit_y,summit_x)*0.5, draw_F(summit_y,summit_x)*0.6, ...
        draw_F(summit_y,summit_x)*0.7, draw_F(summit_y,summit_x)*0.8, ...
        draw_F(summit_y,summit_x)*0.9];%等高線本数
    switch plot_type
        case 'contour'
            contourf(IDSPdata.Vy, IDSPdata.Vx, draw_F, levels)%2次元等高線図
        case 'surf'
            surf(IDSPdata.Vy, IDSPdata.Vx, draw_F)%3次元表面プロット
    end
    colorbar('Ticks',levels,...
        'TickLabels',{'10%','20%','30%','40%','50%','60%','70%','80%','90%'})
    hold on
    for i = 1:n_rows
        plot(IDSPdata.ppoints(i,1), IDSPdata.ppoints(i,2), 'r+')
        hold on
    end
    title(['R = ', num2str(IDSP.r(k,1)),' m'],'FontSize',12);
    daspect([1 1 1])
    xlim([-50,50])
    ylim([-50,50])
    % xticks(-50:10:50)
    % yticks(-50:10:50)
    xlabel('V_{Z} [km/s]','FontSize',12);
    ylabel('V_{R} [km/s]','FontSize',12);
    ax = gca;
    ax.FontSize = 12;
end

if save_fig
    if not(exist([pathname.fig,'/ionvdist'],'dir'))
        mkdir(sprintf("%s", pathname.fig), 'ionvdist');
    end
    saveas(gcf,[pathname.fig,'/ionvdist/','shot', num2str(IDSP.shot),'_',num2str(time),'us_PF1_',num2str(EXP.PF1), ...
        'kV_PF2_',num2str(EXP.PF2),'kV_TF_',num2str(EXP.TF),'kV_EF_',num2str(EXP.EF),'A.png'])
    hold off
    close
else
    hold off
end
end
