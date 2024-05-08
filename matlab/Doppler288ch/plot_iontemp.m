function [] = plot_iontemp(z_IDS,r_IDS,Ti_IDS,date,expval,IDS288ch,pathname,overlay_plot,save_fig)

time = round(IDS288ch.trg+IDS288ch.exp_w/2);%åvë™éûçè

if overlay_plot
    [~,h] = contourf(z_IDS,r_IDS,Ti_IDS,100);
    daspect([1 1 1])
    h.LineStyle = 'none';
    h.FaceAlpha = 0.3;%ìßñæìx(0Ç≈äÆëSÇ…ìßñæ)
    colormap(jet)
    c = colorbar;
    c.Label.String = 'Ion Temperature [eV]';
    if save_fig
        if not(exist([pathname.fig,'/iontemp_psi/',num2str(date)],'dir'))
            mkdir(sprintf("%s/iontemp_psi",pathname.fig), sprintf("%s", num2str(date)));
        end
        saveas(gcf,[pathname.fig,'/iontemp_psi/',num2str(date),'/','shot', num2str(IDS288ch.shot),'_',num2str(time),'us_PF1_',num2str(expval.PF1), ...
            'kV_PF2_',num2str(expval.PF2),'kV_TF_',num2str(expval.TF),'kV_EF_',num2str(expval.EF),'A.png'])
        hold off
        close
    else
        hold off
    end
else
    figure('Position',[600 150 500 400])
    [~,h] = contourf(z_IDS,r_IDS,Ti_IDS,100);
    daspect([1 1 1])
    h.LineStyle = 'none';
    colormap(jet)
    c = colorbar;
    c.Label.String = 'Ion Temperature [eV]';
    title('Local Ion Temperature')
    xlabel('Z [m]')
    ylabel('R [m]')
    if save_fig
        if not(exist([pathname.fig,'/iontemp/',num2str(date)],'dir'))
            mkdir(sprintf("%s/iontemp",pathname.fig), sprintf("%s", num2str(date)));
        end
        saveas(gcf,[pathname.fig,'/iontemp/',num2str(date),'/','shot', num2str(IDS288ch.shot),'_',num2str(time),'us_PF1_',num2str(expval.PF1), ...
            'kV_PF2_',num2str(expval.PF2),'kV_TF_',num2str(expval.TF),'kV_EF_',num2str(expval.EF),'A.png'])
        hold off
        close
    else
        hold off
    end
end

end