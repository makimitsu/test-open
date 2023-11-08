function plot_flow_profile(ExBdata2D,IDSP,FIG,profileplot)
idx_IDSP_z = knnsearch(ExBdata2D.zq(1,:)',IDSP.z(1));
%グラフ
figure('Position', [0 0 500 500],'visible','on');
for i = 1:FIG.tate*FIG.yoko
    switch profileplot
        case 'VExBr'
            plot(ExBdata2D.rq(:,idx_IDSP_z),ExBdata2D.VExB_r(:,idx_IDSP_z,i))
        case 'VExBz'
            plot(ExBdata2D.rq(:,idx_IDSP_z),ExBdata2D.VExB_z(:,idx_IDSP_z,i))
        case '|VExB|'
            plot(ExBdata2D.rq(:,idx_IDSP_z),ExBdata2D.absVExB(:,idx_IDSP_z,i))
    end
    hold on
end
title(sprintf('Z = %0.2f [cm]',ExBdata2D.zq(1,idx_IDSP_z)*1e2))
xlim([0.115 0.275])
xlabel('R [m]')
switch profileplot
    case 'VExBr'
        ylabel('R component of V_{ExB} [km/s]')
    case 'VExBz'
        ylabel('Z component of V_{ExB} [km/s]')
    case '|VExB|'
        ylabel('|V_{ExB}| [km/s]')
end
legendStrings = string(ExBdata2D.time) +"us";
legend(legendStrings)
grid on


end