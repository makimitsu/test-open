function [] = plot_ionvdist(Vx,Vy,F,date,expval,ICCD,pathname,mpoints,ppoints,plot_type,save_fig)
%再構成速度分布関数を描画
[n_rows,~] = size(ppoints);
time = round(ICCD.trg+ICCD.exp_w/2);%プロット時間
for k = 1:mpoints.n_r
    draw_F = reshape(F(:,k), [numel(Vx),numel(Vy)]);%描画用に整形
    draw_F = flipud(draw_F);
    [summit_y,summit_x] = find(draw_F == max(draw_F(:)));
    % fprintf('Summit is (Vz = %.1f, Vr = %.1f)[km/s].\n',Vx(summit_Vx),Vy(summit_Vy));
    levels = [draw_F(summit_y,summit_x)*0.1, draw_F(summit_y,summit_x)*0.2, ...
        draw_F(summit_y,summit_x)*0.3, draw_F(summit_y,summit_x)*0.4, ...
        draw_F(summit_y,summit_x)*0.5, draw_F(summit_y,summit_x)*0.6, ...
        draw_F(summit_y,summit_x)*0.7, draw_F(summit_y,summit_x)*0.8, ...
        draw_F(summit_y,summit_x)*0.9];%等高線本数
    figure('Position',[700 150 550 500],'visible','on')
    switch plot_type
        case 'contour'
            contourf(Vy, Vx, draw_F, levels)%2次元等高線図
        case 'surf'
            surf(Vy, Vx, draw_F)%3次元表面プロット
    end
    colorbar('Ticks',levels,...
        'TickLabels',{'10%','20%','30%','40%','50%','60%','70%','80%','90%'})
    hold on
    for i = 1:n_rows
        plot(ppoints(i,1), ppoints(i,2), 'r+')
        hold on
    end
    % xlim([-50,50])
    % ylim([-50,50])
    xticks(-50:10:50)
    yticks(-50:10:50)
    ax = gca;
    ax.FontSize = 18;
    title(['R = ', num2str(mpoints.r(k,1)),' mm (',num2str(time),' us)'],'FontSize',22);
    xlabel('V_{Z} [km/s]','FontSize',22);
    ylabel('V_{R} [km/s]','FontSize',22);
    % saveas(gcf,[num2str(time),'us_',num2str(r_measured),'mm.png'])
    hold off
    if save_fig
        if not(exist([pathname.fig,'/ionvdist/',num2str(date)],'dir'))
            mkdir(sprintf("%s/ionvdist", pathname.fig), sprintf("%s", num2str(date)));
        end
        saveas(gcf,[pathname.fig,'/ionvdist/',num2str(date),'/','shot', num2str(ICCD.shot),'_',num2str(time),'us_(z,r)_(',num2str(mpoints.z(k,1)),',',num2str(mpoints.r(k,1)),')mm_PF1_',num2str(expval.PF1), ...
            'kV_PF2_',num2str(expval.PF2),'kV_TF_',num2str(expval.TF),'kV_EF_',num2str(expval.EF),'A.png'])
        hold off
        close
    else
        hold off
    end
end
end
